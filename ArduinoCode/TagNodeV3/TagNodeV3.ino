/*
 * Tag:
 * provides ground truth for the localization system
 */

#include <SPI.h>
#include <DW1000.h>
#include <Wire.h>
#include <SparkFunLSM6DS3.h>
#include "Footstep.h"
#define DEVICE_ID       12
#define CLEAR_STEP      true
#define NOT_CLEAR_STEP  false 
#define CLOCK_RATE      48000000
#define SAMPLE_RATE     10

LSM6DS3 myIMU;
uint8_t dataToRead;
int numStep = 1;
int totalStepCount = 0;

// message flow state
volatile byte expectedMsgId = GRANT_TOKEN;
// message sent/received state
volatile boolean sentAck = false;
volatile boolean receivedAck = false;
volatile boolean errorAck = false;
volatile boolean timeoutAck = false;
volatile boolean receiveErrorAck = false;
volatile boolean readDataFlag = false;
volatile boolean alive = false;
volatile byte currentStatus = 0;
volatile byte msgId = 0;
volatile byte msgFrom = 0;
volatile byte msgTo = 0;
volatile byte msgToAnchor = 0;
volatile byte msgToBase = 1;
// data buffer
byte data[LEN_DATA];
byte anchorError[ANCHOR_NUM];
volatile byte accBuffer[LEN_ACC_DATA];
volatile byte bufferIdx = 0;
volatile boolean bufferLooping = false;
// timestamps to remember
DW1000Time timePollSent;
DW1000Time timePollAckReceived;
DW1000Time timeRangeSent;
// watchdog and reset period
unsigned long lastActivity;
unsigned long resetPeriod = 5000;
unsigned long resetSingleEvent = 1000;
unsigned long lastSample;
unsigned long samplePeriod = 100;
// reply times (same on both sides for symm. ranging)
unsigned int replyDelayTimeUS = 3000;
int16_t accX = 0;
int16_t accY = 0;
int16_t accZ = 0;
uint16_t magnitude = 0;
boolean fillBuffer = true;
int selfResetCount = 0;
int timerCheck = 0;
char msg[256];
    

void TC3_Handler(){
  TcCount16* TC = (TcCount16*) TC3; // get timer struct
  if (TC->INTFLAG.bit.OVF == 1) {  // A overflow caused the interrupt
    TC->INTFLAG.bit.OVF = 1;    // writing a one clears the flag ovf flag
    readData();
  }
}

void readData() {
    readDataFlag = true;
}

void setupClock(){
    // Enable clock for TC 
    REG_GCLK_CLKCTRL = (uint16_t) (GCLK_CLKCTRL_CLKEN | GCLK_CLKCTRL_GEN_GCLK0 | GCLK_CLKCTRL_ID_TCC2_TC3) ;
    while ( GCLK->STATUS.bit.SYNCBUSY == 1 ); // wait for sync 
  
    // The type cast must fit with the selected timer mode 
    TcCount16* TC = (TcCount16*) TC3; // get timer struct
  
    TC->CTRLA.reg &= ~TC_CTRLA_ENABLE;   // Disable TC
    while (TC->STATUS.bit.SYNCBUSY == 1); // wait for sync 
  
    TC->CTRLA.reg |= TC_CTRLA_MODE_COUNT16;  // Set Timer counter Mode to 16 bits
    while (TC->STATUS.bit.SYNCBUSY == 1); // wait for sync 
    TC->CTRLA.reg |= TC_CTRLA_WAVEGEN_NFRQ; // Set TC as normal Normal Frq
    while (TC->STATUS.bit.SYNCBUSY == 1); // wait for sync 
  
    TC->CTRLA.reg |= TC_CTRLA_PRESCALER_DIV256;   // Set perscaler
    while (TC->STATUS.bit.SYNCBUSY == 1); // wait for sync 
    
    // TC->PER.reg = 0xFF;   // Set counter Top using the PER register but the 16/32 bit timer counts allway to max  
    // while (TC->STATUS.bit.SYNCBUSY == 1); // wait for sync 
  
    TC->CC[0].reg = 0xFFF;
    while (TC->STATUS.bit.SYNCBUSY == 1); // wait for sync 
    
    // Interrupts 
    TC->INTENSET.reg = 0;              // disable all interrupts
    TC->INTENSET.bit.OVF = 1;          // enable overfollow
    TC->INTENSET.bit.MC0 = 1;          // enable compare match to CC0
  
    // Enable InterruptVector
    NVIC_EnableIRQ(TC3_IRQn);
  
    // Enable TC
    TC->CTRLA.reg |= TC_CTRLA_ENABLE;
    while (TC->STATUS.bit.SYNCBUSY == 1); // wait for sync 

}

void setup() {
    // DEBUG monitoring
    SerialUSB.begin(115200);
    delay(4000);
    SerialUSB.println("### DW1000-arduino-ranging-tag ###");
    // initialize the driver
    DW1000.begin(IRQ_PIN, RESET_PIN);
    DW1000.select(CS_PIN);
    SerialUSB.println("DW1000 initialized ...");
    // general configuration
    DW1000.newConfiguration();
    DW1000.setDefaults();
    DW1000.setDeviceAddress(DEVICE_ID);
    DW1000.setNetworkId(NETWORK_ID);
    DW1000.enableMode(DW1000.MODE_LONGDATA_FAST_ACCURACY);
    DW1000.commitConfiguration();
    SerialUSB.println("Committed configuration ...");
    // DEBUG chip info and registers pretty printed
    DW1000.getPrintableDeviceIdentifier(msg);
    SerialUSB.print("Device ID: "); SerialUSB.println(msg);
    DW1000.getPrintableExtendedUniqueIdentifier(msg);
    SerialUSB.print("Unique ID: "); SerialUSB.println(msg);
    DW1000.getPrintableNetworkIdAndShortAddress(msg);
    SerialUSB.print("Network ID & Device Address: "); SerialUSB.println(msg);
    DW1000.getPrintableDeviceMode(msg);
    SerialUSB.print("Device mode: "); SerialUSB.println(msg);
    // attach callback for (successfully) sent and received messages
    DW1000.attachSentHandler(handleSent);
    DW1000.attachReceivedHandler(handleReceived);
    DW1000.attachErrorHandler(handleError);
    DW1000.attachReceiveFailedHandler(handleReceiveError);
//    DW1000.attachReceiveTimeoutHandler(handleTimeout);
    // anchor starts by transmitting a POLL message
    receiver();
    msgToAnchor = ANCHOR_ID_OFFSET + 1;
    noteActivity();
    
    if( myIMU.begin() != 0 ){
        SerialUSB.println("IMU device error");
    }
    else{
        SerialUSB.println("IMU device OK!");
    }
    //Configure LSM6DS3 as pedometer 
    if( 0 != config_pedometer(CLEAR_STEP) ){
        SerialUSB.println("Configure pedometer fail!");
    }
    SerialUSB.println("Success to Configure pedometer!");
    lastSample = millis();
//    setupClock();
}

void noteActivity() {
    // update activity timestamp, so that we do not reach "resetPeriod"
    lastActivity = millis();
}


void clearErrorFlag(){
    for (int anchorID = 0; anchorID < ANCHOR_NUM; anchorID++){
        anchorError[anchorID] = 0;
    }
}

void resetInactive() {
    // tag sends POLL and listens for POLL_ACK
    currentStatus = TAG_IDLE;
    expectedMsgId = GRANT_TOKEN;
    msgToAnchor = ANCHOR_ID_OFFSET + 1;
    clearDataBuffer();
    clearErrorFlag();
    noteActivity();
}

void handleError(){
    errorAck = true;
}

void handleReceiveError(){
    receiveErrorAck = true;
}

void handleTimeout(){
    timeoutAck = true;
}

void handleSent() {
    // status change on sent success
    sentAck = true;
}

void handleReceived() {
    // status change on received success
    receivedAck = true;
}

void transmitPoll() {
    DW1000.newTransmit();
    DW1000.setDefaults();
    data[0] = POLL;
    data[1] = DEVICE_ID;
    data[2] = msgToAnchor;
    DW1000.setData(data, LEN_DATA);
    DW1000.startTransmit();
}

void transmitTokenAck() {
    SerialUSB.println("send ack");
    DW1000.newTransmit();
    DW1000.setDefaults();
    data[0] = TOKEN_RELEASE;
    data[1] = DEVICE_ID;
    data[2] = msgToBase;
    // TODO: add the vibration data here
    copyAccBuffer();
    DW1000.setData(data, LEN_DATA);
    DW1000.startTransmit();
}

void copyAccBuffer(){
    SerialUSB.println("DEBUG: copy acc to data begin");
    fillBuffer = false;
    if (bufferLooping == false){
        for (int i = 0; i < bufferIdx; i++){
            data[3+i] = accBuffer[i];    
        }  
    } else {
        int dataBufferCounter = 3;
        for (int i = bufferIdx; i < LEN_ACC_DATA; i++){
            data[dataBufferCounter] = accBuffer[i];
            dataBufferCounter += 1;  
        }  
        for (int i = 0; i < bufferIdx; i++){
            data[dataBufferCounter] = accBuffer[i];
            dataBufferCounter += 1;
        }
        bufferLooping = false;
        bufferIdx = 0;
    }
    fillBuffer = true;
    SerialUSB.println("DEBUG: copy acc to data end");
}

void transmitRange() {
    DW1000.newTransmit();
    DW1000.setDefaults();
    data[0] = RANGE;
    data[1] = DEVICE_ID;
    data[2] = msgToAnchor;
    // delay sending the message and remember expected future sent timestamp
    DW1000Time deltaTime = DW1000Time(replyDelayTimeUS, DW_MICROSECONDS);
    timeRangeSent = DW1000.setDelay(deltaTime);
    timePollSent.getTimestamp(data+CONTROL_SIZE);
    timePollAckReceived.getTimestamp(data+CONTROL_SIZE+TSIZE);
    timeRangeSent.getTimestamp(data+CONTROL_SIZE+TSIZE*2);
    DW1000.setData(data, LEN_DATA);
    DW1000.startTransmit();
    //SerialUSB.print("Expect RANGE to be sent @ "); SerialUSB.println(timeRangeSent.getAsFloat());
}

//Setup pedometer mode
int config_pedometer(bool clearStep){
    uint8_t errorAccumulator = 0;
    uint8_t dataToWrite = 0;  //Temporary variable
  
    //Setup the accelerometer******************************
    dataToWrite = 0; 
    
    //  dataToWrite |= LSM6DS3_ACC_GYRO_BW_XL_200Hz;
    dataToWrite |= LSM6DS3_ACC_GYRO_FS_XL_2g;
    dataToWrite |= LSM6DS3_ACC_GYRO_ODR_XL_13Hz;
  
    // Step 1: Configure ODR-26Hz and FS-2g
    errorAccumulator += myIMU.writeRegister(LSM6DS3_ACC_GYRO_CTRL1_XL, dataToWrite);
  
    // Step 2: Set bit Zen_G, Yen_G, Xen_G, FUNC_EN, PEDO_RST_STEP(1 or 0)
    if(clearStep)
      errorAccumulator += myIMU.writeRegister(LSM6DS3_ACC_GYRO_CTRL10_C, 0x3E);
    else
      errorAccumulator += myIMU.writeRegister(LSM6DS3_ACC_GYRO_CTRL10_C, 0x3C);
    
    // Step 3:  Enable pedometer algorithm
    errorAccumulator += myIMU.writeRegister(LSM6DS3_ACC_GYRO_TAP_CFG1, 0x40);
    
    //Step 4: Step Detector interrupt driven to INT1 pin, set bit INT1_FIFO_OVR
    errorAccumulator += myIMU.writeRegister( LSM6DS3_ACC_GYRO_INT1_CTRL, 0x10 );
    
    return errorAccumulator;
}

uint16_t getStepCount(){
    uint8_t dataByte = 0;
    uint16_t stepCount = 0;
    
    myIMU.readRegister(&dataByte, LSM6DS3_ACC_GYRO_STEP_COUNTER_H);
    stepCount = (dataByte << 8) & 0xFFFF;
    
    myIMU.readRegister(&dataByte, LSM6DS3_ACC_GYRO_STEP_COUNTER_L);
    stepCount |=  dataByte;
    return stepCount;
}


void receiver() {
    DW1000.newReceive();
    DW1000.setDefaults();
    // so we don't need to restart the receiver manually
    DW1000.receivePermanently(true);
    DW1000.startReceive();
}

void clearDataBuffer(){
    for (int i = 3; i < LEN_DATA; i++){
        data[i] = 0;  
    }
}

void loop() {
    if (errorAck){
        errorAck = false;
//        DW1000.reset();
        SerialUSB.println("detect error");  
    }

    if (receiveErrorAck){
        receiveErrorAck = false;
        SerialUSB.println("detect receive error");    
    }
    
    if(!sentAck && !receivedAck) {
        // check if inactive
        if(millis() - lastActivity > resetPeriod) {
            SerialUSB.println("DEBUG: self time out");
            resetInactive();
        } else { 
            if (currentStatus == TAG_WAIT_FOR_ANCHOR && millis() - lastActivity > resetSingleEvent){
                SerialUSB.print("Anchor died:");SerialUSB.println(msgToAnchor);
                // mark the tag and continue
                anchorError[msgToAnchor-ANCHOR_ID_OFFSET-1] += 1;
                msgToAnchor += 1;
                if (msgToAnchor > ANCHOR_ID_OFFSET+ANCHOR_NUM){
                    SerialUSB.println("DEBUG: release token");
                    currentStatus = TAG_IDLE;
                    transmitTokenAck();
                    msgToAnchor = ANCHOR_ID_OFFSET+1;
                    expectedMsgId = GRANT_TOKEN;
                } else {
                    delay(ANCHOR_INTERVAL);
                    expectedMsgId = POLL_ACK;
                    transmitPoll();
                }  
                noteActivity();
            }
        }
//        for (int anchorID = 0; anchorID < ANCHOR_NUM; anchorID++){
//            if (anchorError[anchorID] > BS_TAG_RESET_COUNT) {
//                resetInactive();
//            }
//        }
    } else {
        timerCheck = 0;
        selfResetCount = 0;
        // continue on any success confirmation
        if(sentAck) {
            sentAck = false;
            msgId = data[0];
            SerialUSB.print("DEBUG: receive sentAck: ");SerialUSB.println(msgId);
            if(msgId == POLL) {
                currentStatus = TAG_WAIT_FOR_ANCHOR;
                SerialUSB.print("DEBUG: current status: ");SerialUSB.println(currentStatus);
                DW1000.getTransmitTimestamp(timePollSent);
                //SerialUSB.print("Sent POLL @ "); SerialUSB.println(timePollSent.getAsFloat());
            } else if(msgId == RANGE) {
                currentStatus = TAG_WAIT_FOR_ANCHOR;
                DW1000.getTransmitTimestamp(timeRangeSent);
            } else if (msgId == TOKEN_RELEASE) {
                SerialUSB.println("DEBUG: token release sent");
            }
            noteActivity();
            clearDataBuffer();
        }
        if(receivedAck) {
            receivedAck = false;
//            noteActivity();
            // get message and parse
            DW1000.getData(data, LEN_DATA);
            msgId = data[0];
            msgFrom = data[1];
            msgTo = data[2];
            SerialUSB.print("Msg to: "); SerialUSB.print(msgTo);SerialUSB.print(" from: "); SerialUSB.println(msgFrom);
            if (msgTo == DEVICE_ID){
                if(msgId != expectedMsgId) {
                    // unexpected message, start over again
                    SerialUSB.print("ERROR: wrong message: ");SerialUSB.println(msgId);
                    if (currentStatus == TAG_IDLE){
                        // if current status is ranging, ignore this message
                        resetInactive();
                    } else if (msgId != GRANT_TOKEN) {
                        resetInactive();
                    }
                    
                } else {
                    noteActivity();
                    // if it's expected, proceed    
                    if (msgId == GRANT_TOKEN){
                        SerialUSB.println("receive token");
                        clearErrorFlag();
                        expectedMsgId = POLL_ACK;
                        currentStatus = TAG_IDLE;
                        transmitPoll();
                    } else if(msgId == POLL_ACK && msgFrom == msgToAnchor) {
                        SerialUSB.println("DEBUG: receive POLL ACK");
                        DW1000.getReceiveTimestamp(timePollAckReceived);
                        expectedMsgId = RANGE_REPORT;
                        transmitRange();       
                    } else if(msgId == RANGE_REPORT && msgFrom == msgToAnchor) {
                        SerialUSB.println("DEBUG: receive Range Report");
                        currentStatus = TAG_IDLE;
                        float curRange;
                        memcpy(&curRange, data+CONTROL_SIZE, 4);
                        SerialUSB.print(msgFrom);SerialUSB.print(":");SerialUSB.print(curRange);SerialUSB.print("\n");
                        msgToAnchor += 1;
                        if (msgToAnchor > ANCHOR_ID_OFFSET+ANCHOR_NUM){
                            SerialUSB.println("DEBUG: release token");
                            currentStatus = TAG_IDLE;
                            transmitTokenAck();
                            msgToAnchor = ANCHOR_ID_OFFSET+1;
                            expectedMsgId = GRANT_TOKEN;
                        } else {
                            delay(ANCHOR_INTERVAL);
                            expectedMsgId = POLL_ACK;
                            transmitPoll();
                        }
                        noteActivity();
                    } else if(msgId == RANGE_FAILED) {
                        SerialUSB.println("DEBUG: receive ranging failed");
                        resetInactive();
                    }
                    SerialUSB.println("end of processing received");
                }
                
            } else if (msgTo == 0){
                // broadcast
                if(msgId == RESET_NETWORK) {
                    SerialUSB.println("receive request to reset");
                    // unexpected message, start over again
                    resetInactive();
                } 
              
            }
        }
    }

//    if (readDataFlag){
//        readDataFlag = false;
//        accX = myIMU.readRawAccelX();
//        accY = myIMU.readRawAccelY();
//        accZ = myIMU.readRawAccelZ();
//        magnitude = sqrt(accX*accX+accY*accY+accZ*accZ);
//        accBuffer[bufferIdx] = magnitude & 0xFF;
//        bufferIdx += 1;
//        accBuffer[bufferIdx] = (magnitude >> 8) & 0xFF;
//        bufferIdx += 1;
//        if (bufferIdx >= LEN_ACC_DATA){
//            bufferIdx = 0;
//            bufferLooping = true;
//        }
//    }

    if (millis() - lastSample >= samplePeriod){
        accX = myIMU.readRawAccelX();
        accY = myIMU.readRawAccelY();
        accZ = myIMU.readRawAccelZ();
        magnitude = sqrt(accX*accX+accY*accY+accZ*accZ);
        lastSample = millis();
        accBuffer[bufferIdx] = magnitude & 0xFF;
        bufferIdx += 1;
        accBuffer[bufferIdx] = (magnitude >> 8) & 0xFF;
        bufferIdx += 1;
        if (bufferIdx >= LEN_ACC_DATA){
            bufferIdx = 0;
            bufferLooping = true;
        }
    }
    
}

