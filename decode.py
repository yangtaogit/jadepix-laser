#!/usr/bin/env python

'''
MOUDLE : decode source data to cluster with excluding overlap siganl
'''

__author__ = "YANG TAO <yangtao@ihep.ac.cn>"
__copyright__ = "Copyright (c) yangtao"
__created__ = "[2018-04-10 Apr 12:00]"

import os
import sys
import re
import numpy as np
import time
import ROOT
import ConfigParser
import math
import ctypes
# from scipy import ndimage
import logging
logging.basicConfig(level=logging.DEBUG,format= ' %(asctime)s - %(levelname)s- %(message)s')
logging.getLogger().setLevel(logging.INFO)



class Decode():

    def __init__(self,input,output):
        
        self.TEST_FLAG = False

        #define input and output
        self.input = input
        self.output = ROOT.TFile(output,'recreate')

        # configuration
        self.cf = ConfigParser.ConfigParser()
        self.cf.read('configure.conf')

        self.print_number = int(self.cf.get('decode','print_number'))
        self.try_process_number = int(self.cf.get('decode','try_process_number'))
        self.max_frame_number = int(self.cf.get('decode','max_frame_number'))
        self.laser_threshold = int(self.cf.get('decode','laser_threshold'))
        self.binary_threshold = int(self.cf.get('decode','binary_threshold'))

        #init ROOT value

        self.cds_mass_center_list = np.zeros(1, dtype='float32')
        self.binary_mass_center_list = np.zeros(1, dtype='float32')

        self.cds_mass_center_map = ROOT.TH2F('CDS_Mass_Center_Map','CDS_Mass_Center_Map',48,0,48,16,0,16)
        self.binary_mass_center_map =  ROOT.TH2F('Binary_Mass_Center_Map','Binary_Mass_Center_Map',48,0,48,16,0,16)

        self.resolution_tree = ROOT.TTree('Resolution_Tree','Resolution')

        self.resolution_tree.Branch('CDS_Mass_Center_Distribution',self.cds_mass_center_list,'CDS_Mass_Center_Distribution/F')
        self.resolution_tree.Branch('Binary_Mass_Center_Distribution',self.binary_mass_center_list,'Binary_Mass_Center_Distribution/F')


      
    def bytes_to_int(self,thebytes):
        uint8_frame = np.frombuffer(thebytes, dtype=np.uint8)
        uint8_reshape = np.reshape(uint8_frame, newshape=(48, 16, 2))
        logging.debug('\n\n******* reshape uint8 frame *******\n\n')
        logging.debug(str(uint8_reshape))

        int16_frame = np.zeros((48,16), dtype='int16')
        for i in xrange(48):
            for j in xrange(16):
                int16_frame[i][j] = uint8_reshape[i][j][0] + uint8_reshape[i][j][1] * 256
        logging.debug('\n\n******* int16 frame *******\n\n')
        logging.debug(str(int16_frame))

        int_frame = int16_frame.astype('int')
        logging.debug('\n\n******* int frame *******\n\n')
        logging.debug(str(int_frame))

        return int_frame

    def process_raw(self,framebytes):
        clear_frame_bytes = b''
        origin_frame_bytes = framebytes                                          
        tmp = re.findall(b'....(.{32})....',origin_frame_bytes,re.DOTALL)    
        clear_frame_bytes = b''.join(tmp)
        return clear_frame_bytes

    def get_laser_frame(self,cds_frame_adc):
        max_adc = np.max(cds_frame_adc)
        laser_cds_frame = np.zeros((48,16), dtype='int')
        laser_binary_frame = np.zeros((48,16), dtype='int')

        if max_adc > self.laser_threshold :
            laser_cds_frame = cds_frame_adc
            logging.info('\n\n******* signal frame *******\n\n')
            logging.info(str(laser_cds_frame))

            laser_binary_frame = np.where(cds_frame_adc > self.binary_threshold, 1, 0)
            logging.info('\n\n******* binary frame *******\n\n')
            logging.info(str(laser_binary_frame))

            self.get_current_mass_center(laser_cds_frame,laser_binary_frame)

            self.TEST_FLAG = True
            # sys.exit()

        return laser_cds_frame,laser_binary_frame


    def get_current_mass_center(self,laser_cds_frame,laser_binary_frame):
        
        cds_mass_center_up = 0.
        binary_mass_center_up = 0.

        # cds_mass_center = ndimage.measurements.center_of_mass(laser_cds_frame)
        # binary_mass_center = ndimage.measurements.center_of_mass(laser_binary_frame)
        # logging.info('cds center:'+str(cds_mass_center))
        # logging.info('binary center:'+str(binary_mass_center))

        for i in xrange(48):
            for j in xrange(16):
                cds_mass_center_up +=  float(laser_cds_frame[i][j]*math.hypot(i,j))
                binary_mass_center_up +=  float(laser_binary_frame[i][j]*math.hypot(i,j))
        cds_mass_center_down = np.sum(laser_cds_frame)
        binary_mass_center_down = np.sum(laser_binary_frame)

        cds_mass_center = cds_mass_center_up/cds_mass_center_down
        binary_mass_center = binary_mass_center_up/binary_mass_center_down

        self.cds_mass_center_list[0] = cds_mass_center
        self.binary_mass_center_list[0] = binary_mass_center
        self.resolution_tree.Fill()


        logging.info('current cds center:'+str(cds_mass_center))
        logging.info('current binary center:'+str(binary_mass_center))
               
            
    def process_frame(self):

        data =  open(self.input,'rb')
        
        print_number = self.print_number
        try_process_number = self.try_process_number
        max_frame_number = self.try_process_number

        #init value
        seek_position = 0
        frame_number = 0
        cds_frame_number = 0
        broken_frame_number = 0
        broken_bulk_number = 0
        broken_flag = False

        sum_cds_frame = np.zeros((48,16), dtype='int')
        sum_binary_frame = np.zeros((48,16), dtype='int')

        #start 
        while frame_number < max_frame_number:   

            data.seek(seek_position)
            try_process_data = data.read(try_process_number)

            if len(try_process_data) != try_process_number:
                logging.critical('\033[33;1m find total %d frames!\033[0m'%frame_number)
                logging.critical('\033[32;1m find total %d cds frames!\033[0m'%cds_frame_number)
                logging.critical('\033[31;1m find total %d broken frames!\033[0m'%broken_frame_number)
                logging.critical('\033[35;1m find total %d broken bulk!\033[0m'%broken_bulk_number)
                logging.critical(' END !')
                break

            m =re.search(b'(\xaa\xaa\xaa\xaa)(.*?)(\xf0\xf0\xf0\xf0)',try_process_data,re.DOTALL)

            if m:
                if len(m.group(2)) == 1920:
                    frame_number += 1
                    frame_bytes = m.group(2)

                    clear_frame_bytes = self.process_raw(m.group(2))
                    logging.debug('\n\n******** raw frame *******\n\n')
                    logging.debug(repr(clear_frame_bytes))

                    logging.debug('\n\n')

                    frame_adc = self.bytes_to_int(clear_frame_bytes)

                else:
                    data.seek(seek_position+m.start())
                    tmp_process_data = data.read(1928)
                    tmp_m = re.search(b'(\xaa\xaa\xaa\xaa)(.{1920})(\xf0\xf0\xf0\xf0)',tmp_process_data,re.DOTALL)

                    if tmp_m:
                        frame_number += 1
                        frame_bytes = tmp_m.group(2)
                        clear_frame_bytes = self.process_raw(tmp_m.group(2))
                        logging.debug(repr(clear_frame_bytes))

                        logging.debug('\n\n')

                        frame_adc = self.bytes_to_int(clear_frame_bytes)

                    else:
                        broken_frame_number += 1
                        broken_flag = True
                        logging.info('\033[31;1m find %d broken frames!\033[0m'%broken_frame_number)
                        logging.info('\033[31;1m position: (%d %d) \033[0m'%(seek_position+m.start(),seek_position+m.end()))
                        logging.info('\033[31;1m broken length  : %d\033[0m'%len(m.group()))

                ### cds start ####
                if frame_number > 1 :
                    cds_frame_adc = frame_adc-last_frame_adc

                    logging.debug('\n\n******* cds frame *******\n\n')
                    logging.debug(str(cds_frame_adc))
                    laser_cds_frame,laser_binary_frame = self.get_laser_frame(np.negative(cds_frame_adc))
                    
                    cds_frame_number += 1

                    #sum
                    sum_cds_frame += laser_cds_frame
                    sum_binary_frame += laser_binary_frame
                    

                ### cds end ####

                if frame_number % print_number == 0:
                    logging.info('Find %d frames !'%frame_number)
                    logging.info('position: (%d %d)'%(seek_position+m.start(),seek_position+m.end()))
                    logging.info('Get %d cds frames'%cds_frame_number)

                seek_position += (m.start()+(len(m.group())))

                last_frame_adc = frame_adc

            else:
                print('There is no frame in ( %d %d )'%(seek_position,seek_position+try_process_number))
                broken_flag = True
                broken_bulk_number += 1
                seek_position += try_process_number

            if self.TEST_FLAG : 
                 break

        logging.critical('\033[33;1m find total %d frames!\033[0m'%frame_number)
        logging.critical('\033[32;1m find total %d cds frames!\033[0m'%cds_frame_number)
        logging.critical('\033[31;1m find total %d broken frames!\033[0m'%broken_frame_number)
        logging.critical('\033[35;1m find total %d broken bulk!\033[0m'%broken_bulk_number)

        for i in xrange(48):
            for j in xrange(16):
                self.cds_mass_center_map.SetBinContent(i+1,j+1,sum_cds_frame[i][j])
                self.binary_mass_center_map.SetBinContent(i+1,j+1,sum_binary_frame[i][j])

        # self.cds_mass_center_map.Write()
        # self.binary_mass_center_map.Write()
        self.resolution_tree.GetCurrentFile().Write()
        self.output.Close()
        data.close()

    def run(self):
        start_time = time.clock()
        self.process_frame()
        end_time = time.clock()
        print('Running time: %s Seconds'%(end_time-start_time))


