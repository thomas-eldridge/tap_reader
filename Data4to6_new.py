import datetime as dt
import numpy as np
import warnings
import random
import matplotlib.pyplot as plt
from pyorbital.geoloc import ScanGeometry, compute_pixels, get_lonlatalt
import pyorbital.orbital as orb
import pyorbital.astronomy as astro
from find_tle import *
import pdb

# make change to show git diff
# TODO: check data type of every object attribute: are they what you expect?

class Data:
    def __init__(self, the_file):
        """Inputs:
            - the_file; a string representing a path to a .TAP file
        Opens the file to read in binary mode. The file is read and stored in
            - od; the orbit document record (1x per file) containing metadata
              relevant to the whole file
            - dr(s); the data records (arbitrary number per file, usually O(100))
              containing metadata relevant to the upcoming scan block. Each dr contains
              six swath data records (although any number of these may be filled).
        The reader stops when one of several end conditions (in read_header) are met."""
        pointer = open(the_file, 'rb')
        header, end, skip = self.get_header(pointer) # headers are not written on tape, so no endian-ness
        bang = self.zip_bytes_and_goodness(pointer, header)
        self.filename = the_file
        self.od = self.get_orbit_doc(bang, the_file)
        self.dr = []
        footer = self.get_footer(pointer, header)
        i = 0
        while not end:
            header, end, skip = self.get_header(pointer)
            i += 1
            if not end:
                bang = self.zip_bytes_and_goodness(pointer, header)
                footer = self.get_footer(pointer, header)
                if skip:
                    continue
                self.dr.append(self.get_data_rec(bang))
            print i
    def get_header(self, pointer, prev=0):
        """Input:
            - pointer; pointer to the open file
            - prev; indicates the number of zeros encountered already (2x is an EOF mark)
        Output:
            - header; an integer representing the length of the next block to be read
            - EOF; a boolean - True stops the reader, False carries on
            - skip; a boolean - True skips the upcoming data record (including the swath data blocks), False doesn't
        Reads 32 bits from the file as a header. Checks the signbit of the header - if this is on then the upcoming
        block is bad and should be skipped. If the header is zero valued (sometimes happens at the start of the file)
        then another should be read. If two headers are consecutively zero valued then the EOF has been reached. EOF
        can also be reached within a scan block. If this happens, header will read no bytes and return []. This is
        another EOF condition."""
        EOF = False
        skip = False
        raw_header = np.fromfile(pointer, dtype=np.int32, count=1)
        header = raw_header & ((2**31)-1)
        if header != raw_header:
            skip = True
        if header == 0:
            header, EOF, skip = self.get_header(pointer, prev=1)
            if prev == 1:
                EOF = True
        elif len(header) == 0:
            EOF = True
        return header, EOF, skip
    def zip_bytes_and_goodness(self, pointer, header):
        """Inputs:
            - pointer; a pointer to the open TAP file
            - header; an integer representing the number of bytes to be read
        Outputs:
            - bang; Bytes ANd Goodness - a list of the bytes in the upcoming scan block
              zipped together with a boolean representing whether or not each byte passed a check
        A number of bytes is read in, and each byte is put through a parity_check. If the check passes,
        and if the sign bit is off, then the byte is uncorrupted - corrupted otherwise. Returns the list
        of bytes with a second dimension indicating whether or not each byte is OK."""
        head_bytes = np.fromfile(pointer, np.int8, count=header) # reads the whole rest of the file in as bytes
        parity_arr = self.parity(head_bytes) # a companion array for each byte: True if parity is good, else False
        bang = zip(head_bytes, (head_bytes>0)&parity_arr) # zips bytes together with two checks of goodness
        return bang
    def get_footer(self, pointer, header):
        """Inputs:
            - pointer; a pointer to the open TAP file
            - header; an integer representing the number of bytes in a scan block
        Outputs:
            - footer; essentially the same as a header, but referring to the previous scan block
        Uses the get_header method to find the footer. Because headers describe the byte length
        of the upcoming scan block, and footers do so for the previous scan block, each footer should
        match the previous header. If this is not the case the user is warned."""
        footer, end, skip = self.get_header(pointer)
        if header != footer:
            warnings.warn('header and footer do not match')
        return footer
    def parity(self, the_bytes): # could this be sped up?
        """Input:
            - the_bytes; array of 8-bit numbers from the file
        Output:
            - parity_array; array of booleans, corresponding to the_bytes,
              True if the byte passes a parity check, False otherwise
        For each byte, checks to see how many of the least significant six bits are on.
        If the number is even, then the parity bit should be off.
        If the number is odd, then the parity bit should be on.
        If the parity bit is incorrectly assigned, the boolean will be False."""
        parity_array = []
        for byte in the_bytes:
            bits = 0
            good = True
            for bit in range(6):
                if byte & (2**bit) != 0:
                    bits += 1
            if bits % 2 == 0:
                if byte & 0b1000000 == 0:
                    good = False
            else:
                if byte & 0b1000000 != 0:
                    good = False
            parity_array.append(good)
        return parity_array
    def get_orbit_doc(self, bang, the_file):
        """Inputs:
            - bang; the 2D array of bytes (index 0) and goodness (index 1)
            - the_file; a string representing the complete path to the TAP file
        Outputs:
            - Orbit_Doc; an orbit documentation record object for Nimbus 4
        Passes the bytes, their goodness and the filename to the od constructor.
        This method was added because it needs to be overridden in Data2."""
        return Orbit_Doc(bang, the_file)
    def get_data_rec(self, bang):
        """Inputs:
            - bang; the 2D array of bytes (index 0) and goodness (index 1)
        Outputs:
            - Data_Rec; a data record object for Nimbus 4
        Passes the bytes and their goodness to the dr constructor.
        This function was added because it needs to be overridden in Data2."""
        return Data_Rec(bang, self.od)

class Data2(Data):
    def __init__(self, the_file):
        """Inherits from the Data object. Required for Nimbus 5 and 6."""
        Data.__init__(self, the_file)
    def zip_bytes_and_goodness(self, pointer, header):
        """Overrides the method in Data.
        Inputs:
            - pointer; a pointer to the open TAP file
            - header; an integer representing the number of bytes to be read
        Outputs:
            - bang; a list of the bytes in the upcoming scan block
        In the Nimbus 5/6 TAP file, parity and check bits do not exist for each byte. Hence this method
        simply returns the bytes in the upcoming scan block."""
        head_bytes = np.fromfile(pointer, np.int8, count=header) # reads the whole rest of the file in as bytes
        return head_bytes
    def get_orbit_doc(self, bang, the_file):
        """Overrides the method in Data
        Inputs:
            - bang; the array of bytes
            - the_file; a string representing the complete path to the TAP file
        Outputs:
            - Orbit_Doc2; an orbit documentation record object for Nimbus 5 and 6
        Passes the bytes and the filename to the od constructor."""
        return Orbit_Doc2(bang, the_file)
    def get_data_rec(self, bang):
        """Overrides the method in Data
        Inputs:
            - bang; the array of bytes
        Outputs:
            - Data_Rec2; a data record object for Nimbus 5 and 6
        Passes the bytes and their goodness to the dr constructor."""
        return Data_Rec2(bang, self.od)

class Orbit_Doc: # working
    def __init__(self, od_bytes, filename):
        """Inputs:
            - od_bytes; the 2D array of bytes (index 0) and goodness (index 1) for the upcoming scan block
            - filename; a string representing a path to a .TAP file
        Creates an orbit documentation record for the file. The od_bytes are made into 36 bit TAP words in the
        make_words method. Each of these words is then assigned to the relevant attribute of the file.
        Some words must be divided by a scaling factor. The filename is also used in the attribution of the
        file start- and end datetimes, as the bytes passed in have no information on the year of the record."""
        words = self.make_words(od_bytes)
        self.dref = words[0]
        self.nday_start = int(words[2])
        self.start_hour = int(words[3])
        self.start_minute = int(words[4])
        self.start_second = int(words[5])
        self.start_datetime = self.make_start_datetime(filename)
        self.nday_end = int(words[6])
        self.end_hour = int(words[7])
        self.end_minute = int(words[8])
        self.end_second = int(words[9])
        self.end_datetime = self.make_end_datetime(filename)
        if words[10] != -999:
            self.mirror_rot = words[10]/512.
        else:
            self.mirror_rot = words[10] # i.e. -999
        self.sample_freq = words[11]
        self.orbit_no = words[12]
        self.station_code = words[13]
        self.swath_block = words[14]
        self.swaths_per_rec = words[15]
        self.locator_no = words[16]
        self.total_time = self.end_datetime - self.start_datetime
    def make_start_datetime(self, filename):
        """Inputs:
            - filename; the string name of the file
        Makes a datetime object out of the information provided.
        NOTE: This datetime is based on the file name, not the file contents, and should not be trusted over
        the time generated from internal information later in the reader."""
        name = filename.split('/')[-1]
        year = name.split('_')[1][:4]
        int_year = int(year)
        error_in_vars = ((self.nday_start==-999) or (self.start_hour==-999) or (self.start_minute==-999) or (self.start_second==-999))
        if error_in_vars:
            start_datetime = (1970, 01, 01)
        else:
            start_datetime = dt.datetime(int_year, 01, 01)
            start_datetime += dt.timedelta(days=(self.nday_start-1), hours=self.start_hour,
                                      minutes=self.start_minute, seconds=self.start_second)
        return start_datetime
    def make_end_datetime(self, filename):
        """Inputs:
            - filename; the string name of the file
        Makes a datetime object out of the information provided.
        NOTE: This datetime is based on the file name, not the file contents, and should not be trusted over
        the time generated from internal information later in the reader."""
        name = filename.split('/')[-1]
        year = name.split('_')[1][:4]
        int_year = int(year)
        error_in_vars = ((self.nday_end==-999) or (self.end_hour==-999) or (self.end_minute==-999) or (self.end_second==-999))
        if error_in_vars:
            end_datetime = (1970, 01, 01)
        else:
            end_datetime = dt.datetime(int_year, 01, 01)
            end_datetime += dt.timedelta(days=(self.nday_end-1), hours=self.end_hour,
                                      minutes=self.end_minute, seconds=self.end_second)
        return end_datetime
    def make_words(self, the_bytes):
        """Inputs:
            - the_bytes; the 2D array of bytes (index 0) and goodness (index 1)
              for the upcoming scan block
        Outputs:
            - words; a list of TAP words made from the bytes
        Each word is made from 36 bits (4.5 bytes) contained within 48 bits (6 bytes).
        This process is outlined in the docstring for the read_word method."""
        words = []
        for i in range(len(the_bytes)/6):
            words.append(self.read_word(the_bytes[6*i:(6*i)+6]))
        return words
    def read_word(self, word_bytes):
        """Inputs:
            - word_bytes; the six bytes within which a word is buried
        Outputs:
            - word; the word unearthed from the word_bytes
        The THIR recorded data in a six-bit-native-byte machine - hence when read by a modern
        (eight-bit-native-byte) machine two bits are extraneous. These bits are utilised as a check- and parity
        bit, and were handled in the get_bytes_and_goodness method of Data. If the goodness is OK then this method
        proceeds to read the word (otherwise an error value is returned). The six least significant bits from each byte
        are read, and shifted left to make room for the next six. This is occurs a total of six times. The resultant
        word is then returned.
        NOTE: During this process, the endianness of the word is also changed. This is because the THIR machine was
        of opposite endianness to the machine on which this code was written. If this is not the case for the machine
        you are using, this method may need replacing."""
        word = np.int64(0)
        for i in range(len(word_bytes)):
            if not word_bytes[i][1]: # if the goodness is bad for any bit
                word = -999
                break
            else:
                word = word << 6
                mask = word_bytes[i][0] & 0b111111
                word = word | mask
        return word

class Orbit_Doc2(Orbit_Doc): # working
    def __init__(self, od_bytes, filename):
        """Inherits from the Orbit_Doc object. Required for Nimbus 5 and 6."""
        Orbit_Doc.__init__(self, od_bytes, filename)
    def make_words(self, the_bytes):
        """Overrides the method in Orbit_Doc
        Inputs:
            - the_bytes; the array of bytes for the upcoming scan block
        Outputs:
            - words; a list of TAP words made from the bytes
        Each word is made from 36 bits (4.5 bytes). Contrary to the Nimbus 4 TAP files, all bits are
        used for data storage in Nimbus 5 and 6 (hence there is no 'goodness' check). Because 4.5 bytes
        cannot be read, we have no choice but to read 9 bytes and process this as two words back to back.
        This process is outlined in the docstring for the read_word method."""
        words = []
        while len(the_bytes) % 9 != 0:
            the_bytes = np.hstack((the_bytes, [0]))
        for i in range(len(the_bytes)/9):
            word1, word2 = self.read_word(the_bytes[9*i:(9*i)+9])
            words.append(word1)
            words.append(word2)
        return words
    def read_word(self, word_bytes):
        """Overrides the method in Orbit_Doc.
        Inputs:
            - word_bytes; the six bytes within which two words are buried
        Outputs:
            - word1; the word unearthed from the most significant bits in word_bytes
            - word2; the word unearthed from the least significant bits in word_bytes
        See docstring for Orbit_Doc.read_words method. The bytes are added until the word is 40 bits long,
        then shifted back by four bits to produce the first 36 bit word. Then the second 36 bit word is produced
        by picking the 36 least significant bits of the last 40 bits. The resultant words are then returned.
        NOTE: During this process, the endianness of the word is also changed. This is because the THIR machine was
        of opposite endianness to the machine on which this code was written. If this is not the case for the machine
        you are using, this method may need replacing."""
        if len(word_bytes) != 9:
            word1 = -999
            word2 = -999
        else:
            word = np.int64(0)
            for i in range(5):
                word = word << 8
                mask = word_bytes[i] & 0b11111111
                word = word | mask
            word1 = word >> 4
            word = np.int64(0)
            for i in range(4, 9):
                word = word << 8
                mask = word_bytes[i] & 0b11111111
                word = word | mask
            word2 = word & ((2**36)-1)
        return word1, word2

class Data_Rec:
    def __init__(self, dr_bytes, od):
        """Inputs:
            - dr_bytes; the 2D array of bytes (index 0) and goodness (index 1) for the upcoming scan block
            - od; the Orbit_Doc object for this file, containing important metadata
        Creates a data record, containing six swaths. The make_words method returns a list of words to
        use in the definition as well as a marker for when the words should be passed to the Swath_Data constructor.
        Some attributes are set directly (with and without scaling factors), whilst others (namely anchor_nadir_angles
        and sds) are set by methods."""
        words, marker = self.make_words(dr_bytes, od)
        self.nday = words[0]
        self.hour = words[1]
        self.minute = words[2]
        self.second = words[3]
        if words[4] != -999:
            self.roll_error = words[4]/8.
        else:
            self.roll_error = words[4]
        if words[5] != -999:
            self.pitch_error = words[5]/8.
        else:
            self.pitch_error = words[5]
        if words[6] != -999:
            self.yaw_error = words[6]/8.
        else:
            self.yaw_error = words[6]
        self.height = words[7]
        self.cell_temp = words[8]
        self.electro_temp = words[9]
        self.ref_a = words[10]
        self.ref_b = words[11]
        self.ref_c = words[12]
        self.ref_d = words[13]
        self.anchor_nadir_angles = self.set_nadir_angles(words[14:])
        self.sds = self.set_swaths(dr_bytes[marker:], od)
    def make_words(self, the_bytes, od):
        """Inputs:
            - the_bytes; the 2D array of bytes (index 0) and goodness (index 1) for the upcoming scan block
            - od; the Orbit_Doc object for this file, containing important metadata
        Outputs:
            - words; the list of TAP words made from the bytes
            - 6*i+6; the end of the indexes used in the reading of the final word
        Uses the read_half_words and read_full_word methods to produce the list of words. Loops dictate whether
        full_ or half_ words methods are used - this is outlined in the README for the THIR."""
        words = []
        for i in range(7):
            wordA, wordD = self.read_half_words(the_bytes[6*i:(6*i)+6])
            words.append(wordD)
            words.append(wordA)
        for i in range(7, 7+od.locator_no):
            words.append(self.read_full_word(the_bytes[6*i:(6*i)+6]))
        return words, 6*i+6
    def read_half_words(self, some_bytes):
        """Inputs:
            - some_bytes; an array of six bytes to be interpreted
        Outputs:
            - wordA; the least significant 18 bits of a 36 bit word
            - wordD; the most significant 18 bits of a 36 bit word
        Creates a 36 bit word with the read_full_word method, and selects the bits from this to
        make two 18 bit words."""
        word = self.read_full_word(some_bytes)
        if word == -999:
            wordA = -999
            wordD = -999
        else:
            wordA = word & 0x3FFFF
            wordD = word >> 18
        return wordA, wordD
    def read_full_word(self, word_bytes):
        """See docstring for Orbit_Doc.read_word"""
        word = np.int64(0)
        for i in range(len(word_bytes)):
            if not word_bytes[i][1]: # if the goodness is bad for any bit
                word = -999
                break
            else:
                word = word << 6
                mask = word_bytes[i][0] & 0b111111
                word = word | mask
        return word
    def set_nadir_angles(self, words):
        """Inputs:
            - words; a bunch of words chosen to contain the nadir angles
        Outputs:
            - np.array(nadangs)/64.; the resultant nadangs, cast as a numpy array for the sake of easy division
        Loops over the list of words, making sure that signbits register. Returns the nadir angles."""
        nadangs = np.zeros(len(words), dtype=np.int64)
        for i in range(len(words)):
            sign = 1
            if words[i] != -999:
                if (2**35) & words[i] != 0:
                    sign = -1
                nadangs[i] = np.int32(words[i])*sign
                if (nadangs[i] > 4000) | (nadangs[i] < -4000): # out of bounds?
                    nadangs[i] = (-64*999)
            else:
                nadangs[i] = (-64*999) # implicit -999
        return np.array(nadangs)/64.
    def set_swaths(self, sd_bytes, od):
        """Inputs:
            - sd_bytes; the 2D array of bytes (index 0) and goodness (index 1)
            - od; the Orbit_Doc object for this file, containing important metadata
        Outputs:
            - swaths; a list of swath data objects for Nimbus 4
        Passes the bytes and their goodness to the sd constructor six times.
        This function was added because it needs to be overridden in Data_Rec2."""
        swaths = []
        for i in range(od.swaths_per_rec):
            index = i*od.swath_block*6
            swaths.append(Swath_Data(sd_bytes[index:index + (od.swath_block * 6)], od, self))
        return swaths

class Data_Rec2(Data_Rec):
    def __init__(self, dr_bytes, od):
        """Inherits from the Data_Rec object. Required for Nimbus 5 and 6."""
        Data_Rec.__init__(self, dr_bytes, od)
    def make_words(self, the_bytes, od):
        """Overrides the method in Data_Rec
        Inputs:
            - the_bytes; the array of bytes for the upcoming scan block
            - od; the Orbit_Doc object for this file, containing important metadata
        Outputs:
            - words; the list of TAP words made from the bytes
            - 9*i; the end of the indexes used in the reading of the final word
        Uses the read_half_words and read_full_word methods to produce the list of words. Loops dictate whether
        full_ or half_ words methods are used - this is outlined in the README for the THIR."""
        words = []
        for i in range(3):
            word1A, word1D, word2A, word2D = self.read_half_words(the_bytes[9*i:(9*i)+9])
            words.append(word1D)
            words.append(word1A)
            words.append(word2D)
            words.append(word2A)
        # need two more half-words
        word1A, word1D, word2A, word2D = self.read_half_words(the_bytes[27:36])
        words.append(word1D)
        words.append(word1A)
        word1, word2 = self.read_full_word(the_bytes[27:36])
        words.append(word2)
        i = 4
        while len(words) < (14+od.locator_no):
            word1, word2 = self.read_full_word(the_bytes[9*i:(9*i)+9])
            words.append(word1)
            words.append(word2)
            i += 1
        return words, 9*i
    def read_half_words(self, some_bytes):
        """Overrides the method in Data_Rec
        Inputs:
            - some_bytes; an array of nine bytes to be interpreted
        Outputs:
            - word1A; the least significant 18 bits of a 36 bit word
            - word1D; the most significant 18 bits of a 36 bit word
            - word2A; the least significant 18 bits of a 36 bit word
            - word2D; the most significant 18 bits of a 36 bit word
        Creates two 36 bit words with the read_full_word method, and selects the bits from these to
        make four 18 bit words."""
        word1, word2 = self.read_full_word(some_bytes)
        if word1 == -999:
            word1A = -999
            word1D = -999
        else:
            word1A = word1 & 0x3FFFF
            word1D = word1 >> 18
        if word2 == -999:
            word2A = -999
            word2D = -999
        else:
            word2A = word2 & 0x3FFFF
            word2D = word2 >> 18
        return word1A, word1D, word2A, word2D
    def read_full_word(self, word_bytes):
        """Overrides the method in Data_Rec
        See docstring for Orbit_Doc2.read_word"""
        if len(word_bytes) != 9:
            word1 = -999
            word2 = -999
        else:
            word = np.int64(0)
            for i in range(5):
                word = word << 8
                mask = word_bytes[i] & 0b11111111
                word = word | mask
            word1 = word >> 4
            if len(word_bytes) == 9:
                word = np.int64(0)
                mask = word_bytes[i] & 0b1111
                word = word | mask
                for i in range(5,9):
                    word = word << 8
                    mask = word_bytes[i] & 0b11111111
                    word = word | mask
                word2 = word & ((2**36)-1)
            else:
                word2 = -999
        return word1, word2
    def set_swaths(self, sd_bytes, od):
        """Overrides method in Data_Rec
        Inputs:
            - sd_bytes; the array of bytes
            - od; the Orbit_Doc object for this file, containing important metadata
        Outputs:
            - swaths; a list of swath data objects for Nimbus 5 and 6
        Passes the bytes and their goodness to the first sd constructor, then to the second.
        Repeats this process three times. Because of the nature of the Nimbus 5/6 TAP word
        (namely that it is impossible to read a single word) and because the first word of the
        even numbered swaths is always the second half of a couplet, a different constructor is
        required for even numbered swaths."""
        swaths = []
        for i in range(od.swaths_per_rec/2):
            index = i*od.swath_block*9
            swaths.append(Swath_Data2(sd_bytes[index:index+(od.swath_block*9)], od, self))
            swaths.append(Swath_Data3(sd_bytes[index:index+(od.swath_block*9)], od, self))
        return swaths

class Swath_Data:
    def __init__(self, sd_bytes, od, dr):
        """Input:
            - sd_bytes; the bytes relevant for the construction of a swath data record
            - od; the orbit document object for this file, containing important metadata
            - dr; the data record document for this block, containing important metadata
        Creates a swath data object. The make_words method returns a list of words to use in the definition. Some
        of these are set directly (with or without scaling factors), others are set using methods.
        NOTE: If the data population of a given swath is deemed erroneous, that swath will be filled. This should
        obviously be the case when data_pop is zero, but less obviously the data_pop attribute is sometimes corrupted
        and hence enormous. In these cases it is also set to zero and the swath is treated as if it were zero."""
        words = self.make_words(sd_bytes, od)
        if words[0] != -999:
            self.seconds = words[0]/512.
        else:
            self.seconds = words[0]
        self.data_pop = words[1]
        if (self.data_pop != 0) & (self.data_pop < 600):
            if words[2] != -999:
                self.subsat_lat = words[2]/64.
            else:
                self.subsat_lat = words[2]
            if words[3] != -999:
                self.subsat_lon = words[3]/64.
            else:
                self.subsat_lon = words[3]
            self.flags = self.get_flags(words[4])
            self.anchor_lats = self.get_anchor_lats(words[5:5+(2*od.locator_no):2])
            self.anchor_lons = self.get_anchor_lons(words[6:6+(2*od.locator_no):2])
            self.data = self.get_data(words[6+(2*od.locator_no)-1:]) # need to step back by one because of double stepping in lines above
            full_coords = self.get_full_coords(dr, od)
            self.full_lats = full_coords[0]
            self.full_lons = full_coords[1]
        else:
            self.data_pop = 0
            self.subsat_lat = -999
            self.subsat_lon = -999
            array = np.zeros(9)
            array.fill(-999)
            self.flags = array
            array = np.zeros(od.locator_no)
            array.fill(-999)
            self.anchor_lats = np.copy(array)
            self.anchor_lons = np.copy(array)
            self.data = np.array([])
            self.full_lats = np.array([])
            self.full_lons = np.array([])
            self.seconds = -999
    def make_words(self, the_bytes, od):
        """Inputs:
            - the_bytes; the array of bytes for the upcoming swath
            - od; the Orbit_Doc object for this file, containing important metadata
        Outputs:
            - words; the list of TAP words made from the bytes
        Uses the read_half_words and read_full_word methods to produce the list of words. Loops dictate whether
        full_ or half_ words methods are used - this is outlined in the README for the THIR."""
        words = []
        for i in range(2):
            wordA, wordD = self.read_half_words(the_bytes[6*i:(6*i)+6])
            words.append(wordD)
            words.append(wordA)
        j = i+1
        words.append(self.read_full_word(the_bytes[6*j:(6*j)+6]))
        for i in range(j+1, len(the_bytes)):
            wordA, wordD = self.read_half_words(the_bytes[6*i:(6*i)+6])
            words.append(wordD)
            words.append(wordA)
        return words
    def read_half_words(self, some_bytes):
        """See docstring for Data_rec.read_half_words"""
        word = self.read_full_word(some_bytes)
        if word == -999:
            wordA = -999
            wordD = -999
        else:
            wordA = word & 0x3FFFF
            wordD = word >> 18
        return wordA, wordD
    def read_full_word(self, word_bytes):
        """See docstring for Data_rec.read_full_word"""
        word = np.int64(0)
        for i in range(len(word_bytes)):
            if not word_bytes[i][1]: # if the goodness is bad for any bit
                word = -999
                break
            else:
                word = word << 6
                mask = word_bytes[i][0] & 0b111111
                word = word | mask
        return word
    def get_flags(self, number):
        """Inputs:
            - number; an integer to be decoded as flags
        Outputs:
            - flags; the nine flags, allocated an index in an array
        Turns on the relevant flag if the corresponding bit is on in the number.
        Removes unassigned flags before returning. Order: furthest left = flag1."""
        flags = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for i in range(len(flags)):
            if number == -999:
                flags[i] = -999
            else:
                if (number & 2**i) == (2**i):
                    flags[i] = 1
        del flags[12]
        del flags[10]
        del flags[9]
        del flags[6]
        return flags
    def get_anchor_lats(self, latarray):
        """Inputs:
            - latarray; an array of 31 words, corresponding to the latitudes of the anchor points
        Output:
            - retval; the 31 anchor point latitudes, scaled appropriately
        Scales each word in the array."""
        retval = np.zeros(len(latarray))
        for i in range(len(latarray)):
            if (latarray[i] != -999): # raw value error?
                retval[i] = latarray[i]/64.
                if (retval[i] > 180) | (retval[i] < 0): # out of bounds?
                    retval[i] = -999
                else:
                    pass
            else:
                retval[i] = -999
        return retval
    def get_anchor_lons(self, lonarray):
        """Inputs:
            - lonarray; an array of 31 words, corresponding to the longitudes of the anchor points
        Output:
            - retval; the 31 anchor point longitudes, scaled appropriately
        Scales each word in the array."""
        retval = np.zeros(len(lonarray))
        for i in range(len(lonarray)):
            if lonarray[i] != -999:
                retval[i] = lonarray[i]/64.
                if (retval[i] > 360) | (retval[i] < 0): # out of bounds?
                    retval[i] = -999
            else:
                retval[i] = lonarray[i]
        return retval
    def get_data(self, words):
        """Inputs:
            - words; an array of words, corresponding to the data points
        Output:
            - datarray; the data points, scaled appropriately
        Scales each word in the array. The array will have a length dictated by data_pop."""
        if len(words) < self.data_pop:
            datarray = np.array([])
        else:
            datarray = np.zeros(self.data_pop)
            for i in range(self.data_pop):
                if words[i] != -999:
                    datarray[i] = words[i]/8.
                    if (datarray[i] > 400) | (datarray[i] < 0): # out of bounds?
                        datarray[i] = -999
                else:
                    datarray[i] = -999
        return datarray
    def get_full_coords(self, dr, od, tolerance=15):
        """Inputs:
            - dr; the data record document for this block, containing important metadata
            - od; the orbit document object for this file, containing important metadata
            - tolerance; a number describing the extent (in degrees) to which interpolated longitudes may overshoot
              (default=15)
        Outputs:
            - nplats; the interpolated latitudes array
            - nplons; the interpolated longitudes array
        Applies a lagrangian interpolator with a number of points n=5. First, the mirror rotation rate is assumed to be
        constant throughout the flight, and the first (last) anchor point is assumed to correspond to the first
        (last) data point. Then the range of nadir angles is interpolated to provide a value for each data point.
        The 'x' values for the interpolation are the anchor point nadir angles, and the 'y' values are the coordinates.
        The aim is to map our known progression through the x's (the current nadir angle) to the progression through the
        coordinates. For each nadir angle, the find_nearest_index method is used to find the index of the closest angle
        in the anchor points array. Then the surrounding_n method is used to find the n points to use in the
        interpolation. For the sake of interpolating along a smooth function, longitudes are added to/ subtracted from
        if the inspected area crosses the dateline. Then the interp_lagrange method is applied to the relevant info,
        and the output is appended to the array.
        Finally, if longitudes fall within a tolerance of the upper/lower limit then they can be folded back into the
        longitude domain. This compensates for the addition/ subtraction mentioned above. Other erroneous values
        are set as such.
        NOTE: The interpolator works fine for well behaved functions, but the raw anchor point arrays is known to
        behave badly at extreme latitudes. It is this behaviou, not the behaviour of the interpolator, which causes
        a mess in the fully extended coordinate arrays. It is my hope that by fixing the anchor points, the full
        coordinate arrays will behave better and require less 'tolerance'."""
        n = 5
        anchor_nads = dr.anchor_nadir_angles
        full_nads = np.linspace(np.min(anchor_nads), np.max(anchor_nads), self.data_pop)
        lats = []
        lons = []
        for nad in full_nads:
            ind = self.find_nearest_index(anchor_nads, nad)
            relevant_inds = self.surrounding_n(ind, n, od)
            relevant_anchors = anchor_nads[relevant_inds]
            relevant_latitudes = self.anchor_lats[relevant_inds]
            relevant_longitudes = self.anchor_lons[relevant_inds]
            if abs(relevant_longitudes[0]-relevant_longitudes[-1]) > 180:  # i.e. the dateline is crossed
                relevant_longitudes[relevant_longitudes<180]+=360
            lats.append(self.interp_lagrange(nad, relevant_anchors, relevant_latitudes))
            lons.append(self.interp_lagrange(nad, relevant_anchors, relevant_longitudes))
        # deal with slight overshot during interpolation:
        nplons = np.array(lons)
        nplons[(nplons>360)&(nplons<360 + tolerance)] -= 360
        nplons[nplons>360 + tolerance] = -999
        nplons[(nplons>0 - tolerance)&(nplons<0)] += 360
        nplons[nplons<0 - tolerance] = -999
        nplats = np.array(lats)
        nplats[nplats>180] = -999
        nplats[nplats<0] = -999
        return nplats, nplons
    def find_nearest_index(self, the_array, value):
        """Inputs:
            - the_array; the array in which the index is to be found
            - value; the value for which the nearest index in array is to be found
        Output:
            - i; the index where the array is closest to the value
        Subtracts the value across the array and returns the index at which the differences produced are minimised."""
        array = np.copy(the_array)
        array -= value
        i = 0
        while abs(array[i]) != min(abs(array)):
            i += 1
        return i
    def surrounding_n(self, index, n, od):
        """Returns a list of indexes usable for the lagrangian interpolator, depending on the size of n."""
        addresses = np.arange(index-(n/2), index+(n/2)+1)
        while addresses[0] < 0:
            addresses += 1
        while addresses[-1] >= od.locator_no:
            addresses -= 1
        return addresses
    def interp_lagrange(self, x, xs, ys):
        """Inputs:
            - x; the value at which to evaluate the interpolating polynomial
            - xs; the x coordinates of the discrete data points we are trying to interpolate
            - ys; the y coordinates of the discrete data points we are trying to interpolate
        Finds the value of the Lagrange interpolating polynomial through the points
        (xs[i], ys[i]) at x, where 0<=i<=n and n is the length of xs (and therefore of ys).
        Lagrange polynomials which are highly constrained (i.e. n is too large) exhibit
        highly oscillatory behaviour. As such it is recommended that the polynomial be constructed
        based on a few points in the neighborhood."""
        # testing of the function can be done by calling it on a series of points between x_min and x_max and
        # plotting the result alongside a scatter of the original points
        val = 0
        for i in range(len(ys)):
            poly = self.lagrange_poly(x, xs, i)
            val += ys[i]*poly
        return val
    def lagrange_poly(self, x, xs, i):
        """Inputs:
            - x; the x coordinate of the value which interp_lagrange() aims to interpolate
            - xs; the x coordinates of the data points being interpolated by interp_lagrange()
            - i; the index referring to the point in xs about which to create the basis polynomial
            Returns the value of the basis polynomial centred on the i-th coordinate of xs evaluated at x."""
        poly = 1
        for xm in xs[xs!=xs[i]]:
            poly *= (x-xm)/(xs[i]-xm)
        return poly

class Swath_Data2(Swath_Data):
    def __init__(self, sd_bytes, od, dr):
        Swath_Data.__init__(self, sd_bytes, od, dr)
    def make_words(self, the_bytes, od):
        words = []
        word1A, word1D, word2A, word2D = self.read_half_words(the_bytes[:9])
        words.append(word1D)
        words.append(word1A)
        words.append(word2D)
        words.append(word2A)
        word1, word2 = self.read_full_word(the_bytes[9:18])
        words.append(word1)
        word1A, word1D, word2A, word2D = self.read_half_words(the_bytes[9:18])
        words.append(word2D)
        words.append(word2A)
        i = 2
        while len(words) < (5+(2*od.locator_no)+words[1]):
            word1A, word1D, word2A, word2D = self.read_half_words(the_bytes[9*i:(9*i)+9])
            words.append(word1D)
            words.append(word1A)
            words.append(word2D)
            words.append(word2A)
            i += 1
        return words
    def read_half_words(self, some_bytes):
        word1, word2 = self.read_full_word(some_bytes)
        if word1 == -999:
            word1A = -999
            word1D = -999
        else:
            word1A = word1 & 0x3FFFF
            word1D = word1 >> 18
        if word2 == -999:
            word2A = -999
            word2D = -999
        else:
            word2A = word2 & 0x3FFFF
            word2D = word2 >> 18
        return word1A, word1D, word2A, word2D
    def read_full_word(self, word_bytes):
        if len(word_bytes) != 9:
            word1 = -999
            word2 = -999
        else:
            word = np.int64(0)
            for i in range(5):
                word = word << 8
                mask = word_bytes[i] & 0b11111111
                word = word | mask
            word1 = word >> 4
            if len(word_bytes) == 9:
                word = np.int64(0)
                mask = word_bytes[i] & 0b1111
                word = word | mask
                for i in range(4,9):
                    word = word << 8
                    mask = word_bytes[i] & 0b11111111
                    word = word | mask
                word2 = word & ((2**36)-1)
            else:
                word2 = -999
        return word1, word2

class Swath_Data3(Swath_Data2):
    def __init__(self, sd_bytes, od, dr):
        Swath_Data2.__init__(self, sd_bytes, od, dr)
    def make_words(self, the_bytes, od):
        words = []
        midlen = len(the_bytes)/2
        the_other_bytes = np.delete(the_bytes, np.arange(midlen-(9/2)))
        word1A, word1D, word2A, word2D = self.read_half_words(the_other_bytes[:9])
        words.append(word2D)
        words.append(word2A)
        word1A, word1D, word2A, word2D = self.read_half_words(the_other_bytes[9:18])
        words.append(word1D)
        words.append(word1A)
        word1, word2 = self.read_full_word(the_other_bytes[9:18])
        words.append(word2)
        i = 2
        while len(words) < (5+(2*od.locator_no)+words[1]):
            word1A, word1D, word2A, word2D = self.read_half_words(the_other_bytes[9*i:(9*i)+9])
            words.append(word1D)
            words.append(word1A)
            words.append(word2D)
            words.append(word2A)
            i += 1
        return words

class Fields():
    def __init__(self, file_data):
        self.channel = self.get_channel(file_data)
        self.start_time, self.end_time = self.get_time_lims(file_data)
        self.swath_width, self.no_swaths = self.find_swath_dims(file_data)
        self.truetime, self.trueinds, time = self.tdims(file_data) # test this on a more obviously gappy file
        # time is not for recording, but can be used to check that trueinds is working well
        temps = self.set_temps(file_data)
        self.cell_temps = temps[0]
        self.electro_temps = temps[1]
        self.ref_temps_a = temps[2]
        self.ref_temps_b = temps[3]
        self.ref_temps_c = temps[4]
        self.ref_temps_d = temps[5]
        point_errors = self.set_errors(file_data)
        self.roll_errors = point_errors[0]
        self.pitch_errors = point_errors[1]
        self.yaw_errors = point_errors[2]
        the_rest = self.set_rest(file_data)
        self.heights = the_rest[0]
        self.dpops = the_rest[4]
        self.sub_satellite_lats = the_rest[1]
        self.sub_satellite_lons = the_rest[2]
        flags = the_rest[3]
        self.flag1 = flags[:, 0]
        self.flag2 = flags[:, 1]
        self.flag3 = flags[:, 2]
        self.flag4 = flags[:, 3]
        self.flag5 = flags[:, 4]
        self.flag6 = flags[:, 5]
        self.flag8 = flags[:, 6]
        self.flag9 = flags[:, 7]
        self.flag12 = flags[:, 8]
        small_arrays = self.set_small_arrays(file_data)
        self.nadangs = small_arrays[0]
        self.anchor_lats = small_arrays[1]
        self.anchor_lons = small_arrays[2]
        big_arrays = self.set_big_arrays(file_data)
        lons, lats, alts, solzen, solaz, solalt, satzen, sataz = self.geoloc2(file_data)
        self.data = big_arrays[0]
        self.lats = big_arrays[1]
        self.lats2 = lats
        self.lons = big_arrays[2]
        self.lons2 = lons
        self.sol_zen = solzen
        self.sat_zen = satzen
        self.sol_az = solaz
        self.sat_az = sataz
    def find_swath_dims(self, fd):
        widths = []
        for i in range(len(fd.dr)):
            for j in range(len(fd.dr[i].sds)):
                widths.append(fd.dr[i].sds[j].data_pop)
        swaths = len(widths)
        max_width = max(widths)
        return max_width, swaths
    def tdims(self, fd):
        """Finds an array of truetime (len = no_scanlines) which holds the true time for the whole file (dummy times are
        inserted where necessary), as well as an array of trueinds (len = no_obs < no_scanlines) which holds the
        index in truetime of each record in the file."""
        delta0 = self.start_time - dt.datetime(1970,01,01)
        delta1 = self.end_time - dt.datetime(1970,01,01)
        if ((delta0.days*86400) + delta0.seconds > 63072000) & ((delta1.days*86400) + delta1.seconds < 165542400): # then it's a Nimbus 5 file
            delta0 -= dt.timedelta(seconds=12.5)  
            delta1 -= dt.timedelta(seconds=12.5)
        t0 = (delta0.days*86400) + delta0.seconds + (delta0.microseconds/1000000.)
        t1 = (delta1.days*86400) + delta1.seconds + (delta1.microseconds/1000000.)
        scan_sep = 360/fd.od.mirror_rot
        truetime = np.arange(t0, t1+scan_sep, scan_sep)
        trueinds = []
        the_time = []
        for i in range(len(fd.dr)):
            t_base = self.get_tbase(fd, i)
            for j in range(len(fd.dr[i].sds)):
                if fd.dr[i].sds[j].seconds != -999:
                    time = t_base+fd.dr[i].sds[j].seconds-fd.dr[i].second
                    if (time > 63072000) & (time < 165542400): # then Nimbus 5
                        time -= 12.5
                    trueinds.append(np.where(abs(truetime - time)==min(abs(truetime-time)))[0][0])
                    the_time.append(time)
                else:
                    trueinds.append(-1)
                    the_time.append(-999)
        return truetime, trueinds, the_time
    def get_channel(self, obj):
        retval = 'unknown'
        if obj.od.dref == 115:
            retval = 'window'
        elif obj.od.dref == 67:
            retval = 'vapour'
        return retval
    def get_time_lims(self, obj):
        times = []
        for i in range(len(obj.dr)):
            t_base = self.get_tbase(obj, i)
            for j in range(len(obj.dr[i].sds)):
                if obj.dr[i].sds[j].seconds != -999: # times WILL NOT include all entries if there are bad entries!
                    times.append(t_base+obj.dr[i].sds[j].seconds-obj.dr[i].second)
                else:
                    pass
        t0 = min(times)
        time0 = dt.datetime(1970, 01, 01) + dt.timedelta(seconds=t0)
        t1 = max(times)
        time1 = dt.datetime(1970, 01, 01) + dt.timedelta(seconds=t1)
        return time0, time1
    def get_tbase(self, obj, ind):
        day_diff = obj.dr[ind].nday - obj.od.nday_start
        hour_diff = obj.dr[ind].hour - obj.od.start_hour
        min_diff = obj.dr[ind].minute - obj.od.start_minute
        sec_diff = obj.dr[ind].second - obj.od.start_second
        # second is repeated in sd.second, so subtract it from tbase
        seconds = (86400*day_diff) + (3600*hour_diff) + (60*min_diff) + sec_diff
        delta = obj.od.start_datetime - dt.datetime(1970, 01, 01)
        tbase = (delta.days*86400) + (delta.seconds + seconds) + (delta.microseconds/1000000.)
        return tbase
    def set_temps(self, fd):
        array = np.zeros(len(self.truetime))
        array.fill(-999)
        cell_temps = np.copy(array)
        electro_temps = np.copy(array)
        ref_a = np.copy(array)
        ref_b = np.copy(array)
        ref_c = np.copy(array)
        ref_d = np.copy(array)
        for i in range(len(fd.dr)):
            ind = self.trueinds[(fd.od.swaths_per_rec*i)]
            cell_temps[ind] = fd.dr[i].cell_temp
            electro_temps[ind] = fd.dr[i].electro_temp
            ref_a[ind] = fd.dr[i].ref_a
            ref_b[ind] = fd.dr[i].ref_b
            ref_c[ind] = fd.dr[i].ref_c
            ref_d[ind] = fd.dr[i].ref_d
        return (cell_temps, electro_temps, ref_a, ref_b, ref_c, ref_d)
    def set_errors(self, fd):
        array = np.zeros(len(self.truetime))
        array.fill(-999)
        roll_error = np.copy(array)
        pitch_error = np.copy(array)
        yaw_error = np.copy(array)
        for i in range(len(fd.dr)):
            ind = self.trueinds[(fd.od.swaths_per_rec*i)]
            roll_error[ind] = fd.dr[i].roll_error
            pitch_error[ind] = fd.dr[i].pitch_error
            yaw_error[ind] = fd.dr[i].yaw_error
        return (roll_error, pitch_error, yaw_error)
    def set_rest(self, fd):
        array = np.zeros(len(self.truetime))
        array.fill(-999)
        height = np.copy(array)
        sslats = np.copy(array)
        sslons = np.copy(array)
        dpops = np.copy(array)
        flags = np.zeros((len(self.truetime), 9))
        for i in range(len(fd.dr)):
            ind = self.trueinds[(fd.od.swaths_per_rec*i)]
            height[ind] = fd.dr[i].height
            for j in range(len(fd.dr[i].sds)):
                ind = self.trueinds[(fd.od.swaths_per_rec*i)+j]
                sslats[ind] = fd.dr[i].sds[j].subsat_lat
                sslons[ind] = fd.dr[i].sds[j].subsat_lon
                dpops[ind] = fd.dr[i].sds[j].data_pop
                flags[ind] = fd.dr[i].sds[j].flags
        return height, sslats, sslons, flags, dpops
    def set_small_arrays(self, fd):
        array = np.zeros((len(self.truetime), fd.od.locator_no))
        array.fill(-999)
        nads = np.copy(array)
        lats = np.copy(array)
        lons = np.copy(array)
        for i in range(len(fd.dr)):
            for j in range(len(fd.dr[i].sds)):
                ind = self.trueinds[(fd.od.swaths_per_rec*i)+j]
                nads[ind] = fd.dr[i].anchor_nadir_angles
                lats[ind] = fd.dr[i].sds[j].anchor_lats
                lons[ind] = fd.dr[i].sds[j].anchor_lons
        return nads, lats, lons
    def set_big_arrays(self, fd):
        array = np.zeros((len(self.truetime), self.swath_width))
        array.fill(-999)
        data = np.copy(array)
        lats = np.copy(array)
        lons = np.copy(array)
        for i in range(len(fd.dr)):
            for j in range(len(fd.dr[i].sds)):
                ind = self.trueinds[(fd.od.swaths_per_rec*i)+j]
                scanline = self.make_fit(fd.dr[i].sds[j].data)
                latline = self.make_fit(fd.dr[i].sds[j].full_lats)
                lonline = self.make_fit(fd.dr[i].sds[j].full_lons)
                data[ind] = scanline
                lats[ind] = latline
                lons[ind] = lonline
        lats[lats!=-999] -= 90
        lons[lons!=-999] = 360 - lons[lons!=-999]
        lons[(lons!=-999)&(lons>180)] -= 360
        return data, lats, lons
    def make_fit(self, line):
        while len(line) < self.swath_width:
            line = np.hstack((line, [-999]))
        return line
    def geoloc2(self, fd):
        t = self.truetime
        array = np.zeros((len(self.truetime), self.swath_width))
        array.fill(-999)
        lats = np.copy(array)
        lons = np.copy(array)
        alts = np.copy(array)
        sat_zen = np.copy(array)
        sol_zen = np.copy(array)
        sat_az = np.copy(array)
        sol_az = np.copy(array)
        sol_alt = np.copy(array)
        inds = self.trueinds
        roll = self.roll_errors
        pitch = self.pitch_errors
        yaw = self.yaw_errors
        nads = self.nadangs
        dpop = self.dpops
        height = self.heights
        nimbus = fd.filename.split('/')[-1][0] + fd.filename.split('/')[-1][6]
        # find the first non-fill-valued entries for 
        # all relevant components in terms of increasing index:
        current_roll = self.get_initial(roll, inds) - 90
        current_pitch = self.get_initial(pitch, inds) - 90
        current_yaw = self.get_initial(yaw, inds)  - 90
        current_nads = self.get_initial(nads, inds)
        current_pop = self.get_initial(dpop, inds)
        current_height = self.get_initial(height, inds)
        if fd.od.mirror_rot != -999:
            mirror = 360/fd.od.mirror_rot
        else:
            mirror = 1.25
        for i in range(len(t)):
            if roll[i]!=-999:
                current_roll = roll[i] - 90
            if pitch[i]!=-999:
                current_pitch = pitch[i] - 90
            if yaw[i]!=-999:
                current_yaw = yaw[i] - 90
            if all(nads[i]!=-999):
                current_nads = nads[i]
            if (dpop[i]!=-999) & (dpop[i]!=0):
                current_pop = dpop[i]
            if height[i]!=-999:
                current_height = height[i]
            if (i == -1) | (current_pop==0):
                print 'shouldn\'t be here!'
                pos_time = get_geoloc(t[i], current_pop, current_nads, current_roll, current_pitch, current_yaw, nimbus, mirror)
                pos_time.fill(-999)
            else:
                pos_time = get_geoloc(t[i], current_pop, current_nads, current_roll, current_pitch, current_yaw, nimbus, mirror)
            lon = pos_time[0]
            lat = pos_time[1]
            alt = pos_time[2]
            while len(lon) < len(lons[0]):
                lon = np.hstack((lon, -999))
                lat = np.hstack((lat, -999))
                alt = np.hstack((alt, -999))
            lon[np.isnan(lon)] = -999
            lat[np.isnan(lat)] = -999
            alt[np.isnan(alt)] = -999
            lons[i] = lon
            lats[i] = lat
            alts[i] = alt
            print i
            view_angs = np.linspace(current_nads[-1], current_nads[0], current_pop)
            for j in range(len(lon)):
                if (lon[j] == -999) | (lat[j] == -999) | (alt[j] == -999):
                    sol_zen[i, j] = -999
                    sat_zen[i, j] = -999
                    sol_alt[i, j] = -999
                    sol_az[i, j] = -999
                    sat_az[i, j] = -999
                else:
                    now_dt = dt.datetime(1970, 01, 01) + dt.timedelta(seconds=t[i])
                    sol_zen[i, j] = astro.sun_zenith_angle(now_dt, lon[j], lat[j])
                    sol_az[i, j] = np.rad2deg(astro.get_alt_az(now_dt, lon[j], lat[j])[1])
                    # assuming height should be in km
                    the_ind = int(current_pop/2)
                    az, el = orb.get_observer_look(lon[the_ind], lat[the_ind], 
                                                   alt[the_ind], now_dt, lon[j], lat[j], alt[j])
                    sat_zen[i, j] = self.totally_zen(view_angs[j], lat[the_ind], lon[the_ind], lat[j], lon[j])
                    sat_az[i, j] = az
        lons[np.isnan(lons)] = -999
        lats[np.isnan(lats)] = -999
        alts[np.isnan(alts)] = -999
        sat_zen[np.isnan(sat_zen)] = -999
        sol_zen[np.isnan(sol_zen)] = -999
        return lons, lats, alts, sol_zen, sol_az, sol_alt, sat_zen, sat_az
    def get_initial(self, var, inds):
        i = 0
        if type(var[i]) != np.ndarray:
            while var[inds[i]] == -999:
                i += 1
        else:
            while any(var[inds[i]] == -999):
                i += 1
        return var[inds[i]]
    def totally_zen(self, view_ang, sslat, sslon, lat, lon):
        """Inputs:
        - view_ang; the satellite viewing angle to the observed point
        - sslat; the latitude of the subsatellite point
        - sslon; the longitude of the subsatellite point
        - lat; the latitude of the observed point
        - lon; the longitude of the observed point
        Uses the haversine law to calculate the zenith angle, assuming a spherical earth."""
        # convert to radians:
        rsslat = np.deg2rad(sslat)
        rsslon = np.deg2rad(sslon)
        rlat = np.deg2rad(lat)
        rlon = np.deg2rad(lon)
        # find the angle subtended at the earth's centre using the haversine law:
        term1 = np.sin((rlon-rsslon)/2)**2
        term2 = np.cos(rlon)*np.cos(rsslon)*(np.sin((rlat-rsslat)/2)**2)
        theta = 2*np.arcsin(np.sqrt((term1 + term2)))
        dtheta = np.rad2deg(theta)
        # find angle at observation point:
        return (abs(view_ang) + theta)
        
