max=9
for i in `seq 0 $max`
do
echo "$i"
./main "/home/user/Desktop/AT Progs/RandomNetworksFolder/AllNodesNet0$i.txt" "/home/user/Desktop/AT Progs/RandomNetworksFolder/AllScoresNet0$i.txt" "/home/user/Desktop/AT Progs/RandomNetworksFolder/AllConnectionsNet0$i.txt" "/home/user/Desktop/AT Progs/NetTestResults/TestResultNet0$i.txt"
done
