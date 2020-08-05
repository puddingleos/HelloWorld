#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "matplotlibcpp.h"
#include "matprocess.h"


using namespace std;
namespace plt = matplotlibcpp;

bool measureDataRead(string filename,string dataHeadFlag,string dataFlag,unsigned long pointer, vector3DS_t& radarcube){
    ifstream measureData(filename); // read only mode
    if (measureData.is_open()==false)
        return false;

    measureData.seekg(0,ios::end);
    unsigned long fsize = measureData.tellg();// get size of file
    measureData.seekg(0,ios::beg);// move pointer to the beginning

    string headTmp(1025,'\0');
    string dataHead(1025,'\0');
    radarcube.resize(sample_crp, vector2DS_t(crp_frame, vector1DS_t(rx_n,0)));
    //vector3DS_t radarcube3;
    string::size_type indexFlag = 0;
    unsigned long numRead = 0;

    short *dataTemp16 = NULL;
    char *dataTemp = NULL;
    int DATA16;

    while(!measureData.eof()){
        for (string::size_type i=0;i<1024;i++)
            measureData.get(headTmp[i]);
        indexFlag = headTmp.find(dataHeadFlag);
        if (indexFlag!=headTmp.npos){
            int headOffset = indexFlag;
            pointer+=headOffset;
            measureData.seekg(pointer,ios::beg);

            for (string::size_type i=0;i<1024;i++)
                measureData.get(dataHead[i]);
//            pointer+=1024;
            if (dataHead.find(dataFlag)!=dataHead.npos){// find the number of the frame

                for (int i=0;i<30;i++){
                    if (dataHead[i]>='0' && dataHead[i]<'9'){
                        numRead = numRead*10+(dataHead[i]-'0');
                    }
                }
                pointer = measureData.tellg();

                if (numRead == sample_crp*crp_frame*rx_n*2 && pointer+numRead<(unsigned long)fsize){
                    dataTemp = (char*)malloc(sizeof(char)*numRead);
                    measureData.get(dataTemp,numRead);

                    // data_combine
                    dataTemp16 = (short*)malloc(sizeof(short)*numRead/2);
                    for (int ti = 0;ti < numRead/2; ti++){
                        DATA16 = dataTemp[2*ti]*256+dataTemp[2*ti+1];
                        dataTemp16[ti] = (short)(DATA16 > 65536 ? DATA16 - 65536 : DATA16 );
                    }


                    //data_matrix
                    for(int r1=0;r1<sample_crp/8;r1++){
                        for(int r2 = 0; r2 < crp_frame; r2++){
                            for(int r3 = 0; r3 < rx_n; r3++){
                                for (int r4 = 0; r4 < 8; r4++){
                                    radarcube[r1*8+r4][r2][r3] = dataTemp16[r4+8*r3+8*rx_n*r1+rx_n*sample_crp*r2];
                                }
                            }
                        }
                    }


                    //radarcube3.push_back(radarcube);
                    pointer = measureData.tellg();
                    break;
                }
                else{
                    pointer+=numRead;
                    if (pointer<fsize){
                        measureData.seekg(pointer,ios::beg);
                        continue;
                    }
                    else
                        break;
                }
            }
        }
        else{
            if (numRead==0)
                pointer+=1024;
            else
                pointer+=numRead;
            if (pointer<fsize){
                measureData.seekg(pointer,ios::beg);
                continue;
            }
            else
                break;
        }
    }
    measureData.close();// file close
    free(dataTemp);// ptr free
    free(dataTemp16);
//    return radarcube3;
    return true;
}
