/*
 * Copyright (c) 2015 by Thomas Trojer <thomas@trojer.net>
 * Decawave DW1000 library for arduino.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * @file RangingTag.ino
 * Use this to test two-way ranging functionality with two DW1000. This is
 * the tag component's code which polls for range computation. Addressing and
 * frame filtering is currently done in a custom way, as no MAC features are
 * implemented yet.
 *
 * Complements the "RangingAnchor" example sketch.
 */

#include <SPI.h>
#include <DW1000.h>
#include "Footstep.h"
#define DEVICE_ID         1
#define DEBUG             1

// message flow state
//volatile byte expectedMsgId = SYNC_ACK;
volatile byte expectedMsgId = POLL_ACK;
volatile byte currentStatus = BS_IDLE;
// message sent/received state
volatile boolean sentAck = false;
volatile boolean receivedAck = false;
volatile byte msgId = 0;
volatile int msgFrom = 0;
volatile int msgToAnchor = 0;
volatile int msgToTag = 0;

// timestamps to remember
DW1000Time timePollSent;
DW1000Time timePollAckReceived;
DW1000Time timeRangeSent;

// data buffer
byte data[LEN_DATA];
byte accData[LEN_ACC_DATA];
byte tagError[TAG_NUM];
byte anchorError[ANCHOR_NUM];
// watchdog and reset period
unsigned long lastActivity;
unsigned long resetPeriod = 10000;//3000;
unsigned long resetSingleEvent = 2000;
// reply times (same on both sides for symm. ranging)
unsigned int replyDelayTimeUS = 3000;

void noteActivity() {
    // update activity timestamp, so that we do not reach "resetPeriod"
    lastActivity = millis();
}

void resetInactive() {
    // tag sends POLL and listens for POLL_ACK
    SerialUSB.println("reset network");
    networkReset();
    transmitReset();
    clearErrorFlag();
    noteActivity();
}

void clearErrorFlag(){
    for (int tagID = 0; tagID < TAG_NUM; tagID++){
        tagError[tagID] = 0;
    }
    for (int anchorID = 0; anchorID < ANCHOR_NUM; anchorID++){
        anchorError[anchorID] = 0;
    }
}

void networkReset(){
    currentStatus = BS_IDLE;
    msgToAnchor = ANCHOR_ID_OFFSET+1;
    msgToTag = TAG_ID_OFFSET+1;
    expectedMsgId = POLL_ACK;  
//    expectedMsgId = POLL_ACK;            
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

void transmitSync() {
    SerialUSB.println("DEBUG: BROADCAST timestamp");
    DW1000.newTransmit();
    DW1000.setDefaults();
    data[0] = SYNC_REQ;
    data[1] = DEVICE_ID;
    data[2] = 0;
    unsigned long timestamp = micros();
    data[3] = timestamp & 0xFF;
    data[4] = (timestamp>>8) & 0xFF;
    data[5] = (timestamp>>16) & 0xFF;
    data[6] = (timestamp>>24) & 0xFF;
    DW1000.setData(data, LEN_DATA);
    DW1000.startTransmit();
}

void transmitReset() {
    DW1000.newTransmit();
    DW1000.setDefaults();
    data[0] = RESET_NETWORK;
    data[1] = DEVICE_ID;
    data[2] = 0;
    DW1000.setData(data, LEN_DATA);
    DW1000.startTransmit();
}

void transmitToken() {
    DW1000.newTransmit();
    DW1000.setDefaults();
    data[0] = GRANT_TOKEN;
    data[1] = DEVICE_ID;
    data[2] = msgToTag;
    DW1000.setData(data, LEN_DATA);
    DW1000.startTransmit();
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

void receiver() {
    DW1000.newReceive();
    DW1000.setDefaults();
    // so we don't need to restart the receiver manually
    DW1000.receivePermanently(true);
    DW1000.startReceive();
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
    char msg[256];
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
    // anchor starts by transmitting a POLL message
    receiver();
    msgToAnchor = ANCHOR_ID_OFFSET + 1;
    msgToTag = TAG_ID_OFFSET + 1;
//    transmitPoll();
    transmitSync();
    noteActivity();
}

void loop() {
    if (!sentAck && !receivedAck) {
        // check if inactive
        if(millis() - lastActivity > resetPeriod) {
            resetInactive();
        } else {
            // two scenario that will cause reset
            if (currentStatus == BS_WAIT_FOR_TAG && millis() - lastActivity > resetSingleEvent*3) {
                SerialUSB.print("Tag died:");SerialUSB.println(msgToTag);
                // mark the tag and continue
                tagError[msgToTag-TAG_ID_OFFSET-1] += 1;
                msgToTag += 1;
                if (msgToTag > TAG_ID_OFFSET+TAG_NUM){
                    // finished with all tags  
                    delay(BS_INTERVAL);
                    networkReset();
                    transmitSync();
//                    transmitPoll();
                    noteActivity();
                } else {
                    // copy acc data
//                    SerialUSB.write(data+3, LEN_ACC_DATA);
//                    for (int i = 0; i < LEN_ACC_DATA; i++){
//                        accData[i] = data[i+3];  
//                    }
                    expectedMsgId = TOKEN_RELEASE;
                    transmitToken();
                    noteActivity();
                }
            } else if (currentStatus == BS_WAIT_FOR_ANCHOR && millis() - lastActivity >resetSingleEvent){
                SerialUSB.print("Anchor died:");SerialUSB.println(msgToAnchor);
                // mark the tag and continue
                anchorError[msgToAnchor-ANCHOR_ID_OFFSET-1] += 1;
                msgToAnchor += 1;
                if (msgToAnchor > ANCHOR_ID_OFFSET+ANCHOR_NUM){
                    msgToAnchor = ANCHOR_ID_OFFSET+1;
                    // finish anchor round at local
                    // send token to tags       
                    expectedMsgId = TOKEN_RELEASE;
  //                      msgToTag = TAG_ID_OFFSET + 1;
                    transmitToken();
                    noteActivity();
                } else{
                    delay(ANCHOR_INTERVAL);//
//                    expectedMsgId = POLL_ACK;
//                    transmitPoll();
                    expectedMsgId = POLL_ACK;
                    transmitSync();
                    noteActivity();
                }                
            }
        }
        for (int tagID = 0; tagID < TAG_NUM; tagID++){
            if (tagError[tagID] > BS_TAG_RESET_COUNT) {
                resetInactive();
            }
        } 
        
        for (int anchorID = 0; anchorID < ANCHOR_NUM; anchorID++){
            if (anchorError[anchorID] > BS_TAG_RESET_COUNT) {
                resetInactive();
            }
        }
    } else {
        // continue on any success confirmation
        if (sentAck) {
            currentStatus = BS_IDLE;
            sentAck = false;
            msgId = data[0];
            if (msgId == POLL) {
                currentStatus = BS_WAIT_FOR_ANCHOR;
                DW1000.getTransmitTimestamp(timePollSent);
                //SerialUSB.print("Sent POLL @ "); SerialUSB.println(timePollSent.getAsFloat());
            } else if (msgId == RANGE) {
                currentStatus = BS_WAIT_FOR_ANCHOR;
                DW1000.getTransmitTimestamp(timeRangeSent);
                noteActivity();
            } else if (msgId == SYNC_REQ){
                // once finish the sync, start a new poll
                delay(ANCHOR_INTERVAL);
                expectedMsgId = POLL_ACK;
                transmitPoll();
            } else if (msgId == GRANT_TOKEN){
                currentStatus = BS_WAIT_FOR_TAG;
            } else if (msgId == RESET_NETWORK){
                delay(BS_INTERVAL);
                networkReset();
                transmitSync();
            }
            noteActivity();
        }
        if(receivedAck) {
            receivedAck = false;
            // get message and parse
            DW1000.getData(data, LEN_DATA);
            msgId = data[0];
            msgFrom = data[1];
            byte msgTo = data[2];
    //        SerialUSB.print("From:");SerialUSB.print(int(msgFrom));
    //        SerialUSB.print("; To:");SerialUSB.print(int(msgTo));SerialUSB.print("\n");
            if (msgTo == DEVICE_ID){
                currentStatus = BS_IDLE; 
                if(msgId != expectedMsgId) {
                    // unexpected message, start over again
                    SerialUSB.print("error: expect");SerialUSB.print(expectedMsgId);
                    SerialUSB.print("received: ");SerialUSB.println(msgId);
                    networkReset();
                    transmitReset();
                    return;
                } 
                if (msgId == POLL_ACK && msgFrom == msgToAnchor) {
                    DW1000.getReceiveTimestamp(timePollAckReceived);
                    expectedMsgId = RANGE_REPORT;
//                    SerialUSB.println("wait for range report");
                    transmitRange();
                    noteActivity();
                } else if(msgId == RANGE_REPORT && msgFrom == msgToAnchor) {
                    float curRange;
                    memcpy(&curRange, data+CONTROL_SIZE, 4);
                    SerialUSB.print(int(msgFrom));SerialUSB.print(":");SerialUSB.print(curRange);SerialUSB.print("\n");
                    msgToAnchor += 1;
                    if (msgToAnchor > ANCHOR_ID_OFFSET+ANCHOR_NUM){
                      msgToAnchor = ANCHOR_ID_OFFSET+1;
                      // finish anchor round at local
                      // send token to tags       
                      expectedMsgId = TOKEN_RELEASE;
//                      msgToTag = TAG_ID_OFFSET + 1;
                      transmitToken();
                    } else{
                      delay(ANCHOR_INTERVAL);//
                      expectedMsgId = POLL_ACK;
                      transmitPoll();
                    }
                    noteActivity();
                } else if (msgId == TOKEN_RELEASE && msgFrom == msgToTag) {
                    SerialUSB.print("DEBUG: receive token from: "); SerialUSB.println(msgToTag);
                    tagError[msgToTag-TAG_ID_OFFSET-1] = 0;
                    msgToTag += 1;
                    if (msgToTag > TAG_ID_OFFSET+TAG_NUM){
                        // finished with all tags  
                        delay(BS_INTERVAL);
                        networkReset();
//                        transmitPoll();
                        transmitSync();
                        noteActivity();
                    } else {
                        // copy acc data
//                        SerialUSB.write(data+3, LEN_ACC_DATA);
//                        for (int i = 0; i < LEN_ACC_DATA; i++){
//                            accData[i] = data[i+3];  
//                        }
                        delay(TAG_INTERVAL);
                        expectedMsgId = TOKEN_RELEASE;
                        transmitToken();
                        noteActivity();
                    }
                }
                else if(msgId == RANGE_FAILED) {
                    networkReset();
                    transmitPoll();
                    noteActivity();
                }
            }
        }
      
    }
}

