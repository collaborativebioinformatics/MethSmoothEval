#!/usr/bin/env bash

### Remounting the project dirs for updating changes in UI project dir

#Unmount /mnt/project
umount /mnt

#Download dxfuse
wget https://github.com/dnanexus/dxfuse/releases/download/v1.5.0/dxfuse-linux
chmod +x dxfuse-linux

#Create manifest file for dxfuse
echo "{
     \"files\" : [],
     \"directories\" : [
        {
         \"proj_id\" : \"$DX_PROJECT_CONTEXT_ID\",
         \"folder\" : \"/\",
         \"dirname\" : \"/project\"
        }
     ]
   }" > .dxfuse_manifest.json

# Set the parameters:
FUSE_MOUNT="/mnt"
MANIFEST_FILE="/home/dnanexus/.dxfuse_manifest.json"

# Unmount the project:
#sudo umount $FUSE_MOUNT
umount $FUSE_MOUNT

# Restart dxfuse:
/home/dnanexus/dxfuse -readOnly $FUSE_MOUNT $MANIFEST_FILE
#sudo -E /home/dnanexus/dxfuse-linux -readOnly $FUSE_MOUNT $MANIFEST_FILE
