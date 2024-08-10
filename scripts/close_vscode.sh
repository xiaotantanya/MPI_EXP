#!/bin/bash
 
# 获取所有包含vscode关键字的进程ID
pids=$(ps -ef | grep vscode-server| grep -v grep | awk '{print $2}')
 
# 循环遍历并杀死每个进程
for pid in $pids
do
    kill -9 $pid
done