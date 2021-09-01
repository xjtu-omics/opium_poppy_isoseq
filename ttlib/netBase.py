from ctypes import *
import sys
import os
import ftplib

class myFtp:
    ftp = 0

    def __init__(self, host,username,passwd):
        self.ftp = ftplib.FTP(host,username,passwd)

    #download single file
    def downloadFile(self, localFile, remoteFile):
        file_handler = open(localFile,'wb')
        print(file_handler)
        self.ftp.retrbinary(f'RETR {remoteFile}',file_handler.write)
        file_handler.close()
        return True

    #download dir
    def downloadFileTree(self, localDir, remoteDir):
        if not os.path.exists(localDir):
            os.makedirs(localDir)
        self.ftp.cwd(remoteDir)
        remoteNames = self.ftp.nlst()

        for file in remoteNames:
            local = os.path.join(localDir,file)
            try:
                self.downloadFile(local,file)
            except:
                if not os.path.isdir(local):
                    os.remove(local)
                    os.makedirs(local)
                self.downloadFileTree(local,file)
        self.ftp.cwd("..")
        return

    def close(self):
        self.ftp.quit()