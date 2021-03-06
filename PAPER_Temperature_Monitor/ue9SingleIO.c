//Author: LabJack
//Oct. 2, 2007
//This example program makes 3 SingleIO low-level function calls.  One call sets
//DAC0 to 2.500 V.  One call reads voltage from AIN0.  One call reads the 
//temperature from the internal temperature sensor.  Control firmware version
//1.03 and above needed for SingleIO.
#include "ue9.h"

const int ue9_port = 52360;

int singleIO_AV_example(int socketFD, ue9CalibrationInfo *caliInfo);

int main(int argc, char **argv)
{
  int socketFD;
  ue9CalibrationInfo caliInfo;
  int loop;
  //long start, total;

  if(argc < 2)
  {
    printf("Please enter an ip address to connect to.\n");
    exit(0);
  }
  else if(argc > 2)
  {
    printf("Too many arguments.\nPlease enter only an ip address.\n");
    exit(0);
  }

  //Opening TCP connection to UE9
  if( (socketFD = openTCPConnection(argv[1], ue9_port)) < 0)
    goto done;

  //Getting calibration information from UE9
  if(getCalibrationInfo(socketFD, &caliInfo) < 0)
    goto close;

  //start = getTickCount();
  for(loop = 0; loop < 1000; loop++)
    {
  if(singleIO_AV_example(socketFD, &caliInfo) < 0)
    goto close;
    }

  //total = getTickCount() - start;

  //printf("1000 IOStream reads in %ld ms.\n", total);

close:
  if(closeTCPConnection(socketFD) < 0)
    printf("Error: failed to close socket\n");
done:
  return 0;
}

//Sends 3 SingleIO low-level commands to set DAC0, read AIN0 and read the 
//temperature
int singleIO_AV_example(int socketFD, ue9CalibrationInfo *caliInfo)
{
  uint8 sendBuff[8], recBuff[8];
  int sendChars, recChars;
  double voltage;
  double temperature; //in Kelvins
  uint16 bytesVoltage, bytesTemperature;
  uint8 ainResolution;

  ainResolution = 12;
  //ainResolution = 18;  //high-res mode for UE9 Pro only

  /* set voltage of DAC0 to 2.500 V */
  //if(analogToUncalibratedBinaryVoltage(2.500, &bytesVoltage) < 0)  
  if(analogToCalibratedBinaryVoltage(caliInfo, 0, 2.500, &bytesVoltage) < 0) 
    return -1;

  sendBuff[1] = (uint8)(0xA3);  //command byte
  sendBuff[2] = (uint8)(0x05);  //IOType = 5 (analog)
  sendBuff[3] = (uint8)(0x00);  //Channel = 0 (DAC0)
  sendBuff[4] = (uint8)( bytesVoltage & (0x00FF) );  //low bits of voltage 
  sendBuff[5] = (uint8)(bytesVoltage / 256) + 192;    //high bits of voltage
                                                     //(bit 7 : Enable, bit 6: Update)
  sendBuff[6] = (uint8)(0x00);  //Settling time - does not apply to analog output
  sendBuff[7] = (uint8)(0x00);  //reserved
  sendBuff[0] = normalChecksum8(sendBuff, 8);

  //Sending command to UE9
  sendChars = send(socketFD, sendBuff, 8, 0);
  if(sendChars < 8)
  {
    if(sendChars == -1)
      goto sendError0;
    else
      goto sendError1;
  }

  //Receiving response from UE9
  recChars = recv(socketFD, recBuff, 8, 0);
  if(recChars < 8)
  {
    if(recChars == -1)
      goto recvError0;
    else
      goto recvError1;
  }

  if((uint8)(normalChecksum8(recBuff, 8)) != recBuff[0])
    goto chksumError;

  if(recBuff[1] != (uint8)(0xA3))
    goto commandByteError;

  /*
  if(recBuff[2] != (uint8)(0x05))
    goto IOTypeError;

  if(recBuff[3] != 0)
    goto channelError;
  */

  printf("Set DAC0 voltage to 2.500 V ...\n");

  /* read voltage from AIN0 */
  sendBuff[1] = (uint8)(0xA3);  //command byte
  sendBuff[2] = (uint8)(0x04);  //IOType = 4
  sendBuff[3] = (uint8)(0x00);  //Channel = 0 (AIN0)
  sendBuff[4] = (uint8)(0x00);  //BipGain (Bip = unipolar, Gain = 1)
  sendBuff[5] = ainResolution;   //Resolution = 12
  sendBuff[6] = (uint8)(0x00);  //SettlingTime = 0
  sendBuff[7] = (uint8)(0x00);  //Reserved
  sendBuff[0] = normalChecksum8(sendBuff, 8);

  //Sending command to UE9
  sendChars = send(socketFD, sendBuff, 8, 0);
  if(sendChars < 8)
  {
    if(sendChars == -1)
      goto sendError0;
    else  
      goto sendError1;
  }

  //Receiving response from UE9
  recChars = recv(socketFD, recBuff, 8, 0);
  if(recChars < 8)
  {
    if(recChars == -1)
      goto recvError0;
    else  
      goto recvError1;
  }

  if((uint8)(normalChecksum8(recBuff, 8)) != recBuff[0])
    goto chksumError;

  if(recBuff[1] != (uint8)(0xA3))
    goto commandByteError;

  if(recBuff[2] != (uint8)(0x04))
    goto IOTypeError;

  if(recBuff[3] != 0)
    goto channelError;

  bytesVoltage = recBuff[5] + recBuff[6] * 256;

  //if(binaryToUncalibratedAnalogVoltage(sendBuff[4], bytesVoltage, &voltage) < 0)
  if(binaryToCalibratedAnalogVoltage(caliInfo, sendBuff[4], ainResolution, bytesVoltage, &voltage) < 0) 
    return -1;

  printf("Voltage read from AI0: %.6f V\n", voltage);

  /* read temperature from internal temperature sensor */
  sendBuff[1] = (uint8)(0xA3);  //command byte
  sendBuff[2] = (uint8)(0x04);  //IOType = 4 (analog in)
  sendBuff[3] = (uint8)(0x85);  //Channel = 133 (tempSensor)
  sendBuff[4] = (uint8)(0x00);  //Gain = 1 (Bip does not apply)  
  sendBuff[5] = (uint8)(0x0C);  //Resolution = 12
  sendBuff[6] = (uint8)(0x00);  //SettlingTime = 0
  sendBuff[7] = (uint8)(0x00);  //Reserved
  sendBuff[0] = normalChecksum8(sendBuff, 8);

  //Sending command to UE9
  sendChars = send(socketFD, sendBuff, 8, 0);
  if(sendChars < 8)
  {
    if(sendChars == -1)
      goto sendError0;
    else  
      goto sendError1;
  }

  //Receiving response from UE9
  recChars = recv(socketFD, recBuff, 8, 0);
  if(recChars < 8)
  {
    if(recChars == -1)
      goto recvError0;
    else  
      goto recvError1;
  }

  if((uint8)(normalChecksum8(recBuff, 8)) != recBuff[0])
    goto chksumError;

  if(recBuff[1] != (uint8)(0xA3))
    goto commandByteError;

  if(recBuff[2] != (uint8)(0x04))
    goto IOTypeError;

  if(recBuff[3] != (uint8)(0x85))
    goto channelError;

  bytesTemperature = recBuff[5] + recBuff[6] * 256;

  //assuming high power level
  if(binaryToCalibratedAnalogTemperature(caliInfo, 0, bytesTemperature, &temperature) < 0)
    return -1;

  printf("Temperature read internal temperature sensor (channel 133): %.1f K\n\n", temperature);
  return 0;

//error printouts
sendError0:
  printf("Error : send failed\n");
  return -1;
sendError1:
  printf("Error : did not send all of the buffer\n");
  return -1;
recvError0:
  printf("Error : recv failed\n");
  return -1;
recvError1:  
  printf("Error : recv did not receive all of the buffer\n");
  return -1;
chksumError:
  printf("Error : received buffer has bad checksum\n");
  return -1;
commandByteError:
  printf("Error : received buffer has wrong command byte\n");
  return -1;
IOTypeError:  
  printf("Error : received buffer has wrong IOType\n");
  return -1;
channelError:  
  printf("Error : received buffer has wrong channel\n");
  return -1;
}
