Universe = vanilla
Executable = monopoles.sh
Initialdir = /net/students/home/as4ae/monopoles
#Log = Monopoles.$(Cluster).$(Process).log
#Output = Monopoles.$(Cluster).$(Process).txt
#Error = Monopoles.$(Cluster).$(Process).error
getenv = TRUE
arguments = $(Cluster) $(Process) /net/students/home/as4ae/monopoles
notification = NEVER
queue 10
