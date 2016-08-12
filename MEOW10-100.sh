max=99
for i in `seq 10 $max`
do
echo "$i"
./main "/home/user/Desktop/AT Progs/RandomNetworksFolder/AllNodesNet$i.txt" "/home/user/Desktop/AT Progs/RandomNetworksFolder/AllScoresNet$i.txt" "/home/user/Desktop/AT Progs/RandomNetworksFolder/AllConnectionsNet$i.txt" "/home/user/Desktop/AT Progs/NetTestResults/TestResultNet$i.txt"
done
