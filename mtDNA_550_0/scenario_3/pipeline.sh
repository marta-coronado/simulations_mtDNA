#!/bin/bash 
#title          :pipeline.sh
#description    :This script perform simulations with Recodon v1.6.0 from a configuration 
                #file and execute PhiPack to detect possible events of recRate.
#author         :Marta Coronado Zamora
#date           :2018-12-02
#version        :0.1      
#usage          :./pipeline.sh     
#bash_version   :4.2.46(2)-release
#============================================================================


### DECLARE VARIABLES ###

# Parameters from config recodon
scenarioType='' # parameter combination defining an scenario
sequences='' # number of sequences to simulate
popSize='' # effective population size
recRate='' # recombination rate
seed='' # replicability

# Invariable parameters
numReplicates=100 # number of replicates
length=16478 # sequence length
mutRate=4.1e-7 # mutation rate
invSites=0 # invariable sites
status=1 # diploid (2) or haploid (1)

# Parameters defining the scenario type
popType=0 #constant population
popGrowth=0 #population growth
alphaShape=0 #gamma distribution for heterogeneous mutation rate

# Jukes and Cantor model
nucPropoportion="0.25 0.25 0.25 0.25" #proportion of nucleotides
transProportion=0.5 #transversion and transition probability 
jukesCantor=0

if [[ $nucPropoportion == "0.25 0.25 0.25 0.25" && $transProportion == "0.5" ]]; then
	jukesCantor=1
fi

# Output
file_type=2 #fasta
info=1 #level of data displayed in the console

mkdir Simulations
> Simulations/PhiPack_global_result.txt

### READ CONFIG_RECODON FILE

while IFS='' read -r scenario || [[ -n "$scenario" ]]; do

    scenarioType=$(echo $scenario | cut -f1 -d' ')
    sequences=$(echo $scenario | cut -f2 -d' ')
    popSize=$(echo $scenario | cut -f3 -d' ')
    recRate=$(echo $scenario | cut -f4 -d' ')
    seed=$(echo $scenario | cut -f5 -d' ')

    if [[ $scenarioType == 1 ]]; then
    	echo "Scenario 1"
    	
    	./Recodon-1.6.0/exe/recodon1.6.0 -n$numReplicates -s$sequences -l$length -e$popSize -p$popType -r$recRate -u$mutRate -f4 $nucPropoportion -t$transProportion -i$invSites -_$status -y$info -*$file_type -#$seed
    
    	mkdir Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate

   		mv Results/sequences Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta
   		
   		sequences_split=$(echo "($sequences*2)+1" | bc)
   		split -l $sequences_split -a 3 --numeric-suffixes=1 --additional-suffix=.fasta  Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r${recRate}_
   		
   		mkdir Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/FASTA
   		mv Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/FASTA/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta
   		
   		> Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt

   		for simulations in `ls Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario*`
   			do
   				simulation=$(echo $simulations | cut -f3 -d'/' | cut -f1,2,3,4 -d'_')
   				numSimulation=$(echo $simulations | cut -f3 -d'/' | cut -f5 -d'_' | cut -f1 -d'.' | sed 's/^0*//')

   				./PhiPack/Phi -o -f $simulations
   				pValNSS=$(grep "NSS:" Phi.log | tr -s ' ' | cut -f2 -d' ') 
   				pValMaxChi=$(grep "Max Chi^2:" Phi.log | tr -s ' ' | cut -f3 -d' ') 
   				pValPHI=$(grep "PHI (Normal):" Phi.log | tr -s ' ' | cut -f3 -d' ') 

   				rho=$(echo | awk "BEGIN { print $popSize * 2 * $recRate * $length }")
   				thetaRate=$(echo | awk "BEGIN { print $popSize * $mutRate }")
   				theta=$(echo | awk "BEGIN { print $popSize * $mutRate * $length }")
		
   				echo $scenarioType $simulation $numSimulation $sequences $popSize $recRate $mutRate $length $rho $thetaRate $theta $status $invSites $jukesCantor $popGrowth $alphaShape $pValNSS $pValMaxChi $pValPHI >> Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt
   				echo $scenarioType $simulation $numSimulation $sequences $popSize $recRate $mutRate $length $rho $thetaRate $theta $status $invSites $jukesCantor $popGrowth $alphaShape $pValNSS $pValMaxChi $pValPHI >> Simulations/PhiPack_global_result.txt

   			done

		sed  -i '1i scenarioType simulation numSimulation numSequences popSize recRate mutRate length rho thetaRate theta status invSites jukesCantor popGrowth alphaShape pValNSS pValMaxChi pValPHI' Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt

    fi

    if [[ $scenarioType == 2 ]]; then
    	echo "Scenario 2"
    	alphaShape=0.2
    	./Recodon-1.6.0/exe/recodon1.6.0 -n$numReplicates -s$sequences -l$length -e$popSize -p$popType -r$recRate -u$mutRate -a$alphaShape -f4 $nucPropoportion -t$transProportion -i$invSites -_$status -y$info -*$file_type -#$seed
    
    	mkdir Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate

   		mv Results/sequences Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta
   		
   		sequences_split=$(echo "($sequences*2)+1" | bc)
   		split -l $sequences_split -a 3 --numeric-suffixes=1 --additional-suffix=.fasta  Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r${recRate}_

   		mkdir Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/FASTA
   		mv Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/FASTA/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta

  		> Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt

   		for simulations in `ls Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario*`
   			do
 
 				simulation=$(echo $simulations | cut -f3 -d'/' | cut -f1,2,3,4 -d'_')
   				numSimulation=$(echo $simulations | cut -f3 -d'/' | cut -f5 -d'_' | cut -f1 -d'.' | sed 's/^0*//')

   				./PhiPack/Phi -o -f $simulations
   				pValNSS=$(grep "NSS:" Phi.log | tr -s ' ' | cut -f2 -d' ') 
   				pValMaxChi=$(grep "Max Chi^2:" Phi.log | tr -s ' ' | cut -f3 -d' ') 
   				pValPHI=$(grep "PHI (Normal):" Phi.log | tr -s ' ' | cut -f3 -d' ') 

   				rho=$(echo | awk "BEGIN { print $popSize * 2 * $recRate * $length }")
   				thetaRate=$(echo | awk "BEGIN { print $popSize * $mutRate }")
   				theta=$(echo | awk "BEGIN { print $popSize * $mutRate * $length }")
			
   				echo $scenarioType $simulation $numSimulation $sequences $popSize $recRate $mutRate $length $rho $thetaRate $theta $status $invSites $jukesCantor $popGrowth $alphaShape $pValNSS $pValMaxChi $pValPHI >> Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt
   				echo $scenarioType $simulation $numSimulation $sequences $popSize $recRate $mutRate $length $rho $thetaRate $theta $status $invSites $jukesCantor $popGrowth $alphaShape $pValNSS $pValMaxChi $pValPHI >> Simulations/PhiPack_global_result.txt
   			done

		sed  -i '1i scenarioType simulation numSimulation numSequences popSize recRate mutRate length rho thetaRate theta status invSites jukesCantor popGrowth alphaShape pValNSS pValMaxChi pValPHI' Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt

    fi

    if [[ $scenarioType == 3 ]]; then
    	echo "Scenario 3"
    	
    	popGrowth=1e-3
    	alphaShape=0.2

    	./Recodon-1.6.0/exe/recodon1.6.0 -n$numReplicates -s$sequences -l$length -e$popSize -g$popGrowth -r$recRate -u$mutRate -a$alphaShape -f4 $nucPropoportion -t$transProportion -i$invSites -_$status -y$info -*$file_type -#$seed
    
    	mkdir Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate

   		mv Results/sequences Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta
   		
   		sequences_split=$(echo "($sequences*2)+1" | bc)
   		split -l $sequences_split -a 3 --numeric-suffixes=1 --additional-suffix=.fasta  Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r${recRate}_

   		mkdir Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/FASTA
   		mv Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/FASTA/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate.fasta
   		
  		> Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt

   		for simulations in `ls Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/scenario*`
   			do
   				simulation=$(echo $simulations | cut -f3 -d'/' | cut -f1,2,3,4 -d'_')
   				numSimulation=$(echo $simulations | cut -f3 -d'/' | cut -f5 -d'_' | cut -f1 -d'.' | sed 's/^0*//' )

   				./PhiPack/Phi -o -f $simulations
   				pValNSS=$(grep "NSS:" Phi.log | tr -s ' ' | cut -f2 -d' ') 
   				pValMaxChi=$(grep "Max Chi^2:" Phi.log | tr -s ' ' | cut -f3 -d' ') 
   				pValPHI=$(grep "PHI (Normal):" Phi.log | tr -s ' ' | cut -f3 -d' ') 

   				rho=$(echo | awk "BEGIN { print $popSize * 2 * $recRate * $length }")
   				thetaRate=$(echo | awk "BEGIN { print $popSize * $mutRate }")
   				theta=$(echo | awk "BEGIN { print $popSize * $mutRate * $length }")

   				
   				echo $scenarioType $simulation $numSimulation $sequences $popSize $recRate $mutRate $length $rho $thetaRate $theta $status $invSites $jukesCantor $popGrowth $alphaShape $pValNSS $pValMaxChi $pValPHI >> Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt
   				echo $scenarioType $simulation $numSimulation $sequences $popSize $recRate $mutRate $length $rho $thetaRate $theta $status $invSites $jukesCantor $popGrowth $alphaShape $pValNSS $pValMaxChi $pValPHI >> Simulations/PhiPack_global_result.txt
   			done
   		sed  -i '1i scenarioType simulation numSimulation numSequences popSize recRate mutRate length rho thetaRate theta status invSites jukesCantor popGrowth alphaShape pValNSS pValMaxChi pValPHI' Simulations/scenario${scenarioType}_n${sequences}_popSize${popSize}_r$recRate/PhiPack_result.txt
    fi

sequences_split=0

done < <(tail -n "+2" config_recodon) # read the config file skipping the header

sed  -i '1i scenarioType simulation numSimulation numSequences popSize recRate mutRate length rho thetaRate theta status invSites jukesCantor popGrowth alphaShape pValNSS pValMaxChi pValPHI' Simulations/PhiPack_global_result.txt

rm Phi.inf.list
rm Phi.inf.sites
rm Phi.log
rm Phi.poly.unambig.sites
rm -r Results