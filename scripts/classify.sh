# Examples:
# SVM train classification:
# ./classify.sh ../datasets/Synthetic/synthetic.sparse 1 1 ../ ../ "" "-n 4" 2 0 0
# SVM test classification:
# ./classify.sh ../datasets/Synthetic/synthetic.sparse 2 1 ../ ../ "" "-n 4" 2 0 0
# Random Forest train then test
# ./classify.sh ../datasets/CLData/test/dev-first10_.arff 3 0 ../datasets/CLData/train/train-first10_.arff ./RF-out.out "" "-n 4" 3 10 10

# apply parallel classification on data
# We support SVM ad Random Forest
# script parameters are:
# 1- input file
# 2- operation type; train or classify: 1-SVM train, 2-SVM classify, 3- RF train and classify
# 3- used clasification model
#	1- SVM
#	2- Random Forest
# 4- model (path for SVM, training file for RF)
# 5- output file
# 6- application specific parameters
# 7- MPI specific parameters
# 8- file type: 0-csv, 1-csv with header, 2-  sparse csv, 3- ARFF file, 4- binary (4-bytes), 5- binary (8-bytes)
# 9- classification parameter index [from]
# 10- classification parameter index [to]

inputFile=$1

# Train
if [ $2 -eq 1 ]
then
	echo "SVM Training ..."
	case $3 in
		1)
		    # SVM
		    if [ $8 -ne 2 ]
                    then
			echo "File needs conversion ..."
			inputFile=tempooo
			../systems/FileConverter/converter $1 $inputFile 0 $8 2 $9 ${10}
		    fi
		    mpiexec $7 ../systems/ParallelSVM/psvm/svm_train -model_path $4 $inputFile
		    ;;
		*)
		    echo "Usage: parameter #2 should be 1"
		    exit 1
		    
	esac
fi
if [ $2 -eq 2 ]
then
	echo "SVM Testing ..."
	case $3 in
                1)
                    # SVM
		    if [ $8 -ne 2 ]
                    then
                        echo "File needs conversion ..."
                        inputFile=tempooo
                        ../systems/FileConverter/converter $1 $inputFile 0 $8 2 $9 ${10}
                    fi
                    mpiexec $7 ../systems/ParallelSVM/psvm/svm_predict -model_path $4 -output_path $5 $inputFile
		    echo "Output file is: "$5"PredictResult"
                    ;;
                *)
                    echo "Usage: parameter #2 should be between 1-2"
                    exit 1

	esac
fi
if [ $2 -eq 3 ]
then
	trainFile=$4
	echo "Random Forest training and testing"
	if [ $8 -ne 2 ]
        then
        	echo "Files needs conversion ..."
		echo "Converting testing file ..."
                inputFile=tempooo
                ../systems/FileConverter/converter $1 $inputFile 0 $8 3 $9 ${10}
		echo "converting training file ..."
		trainFile=tempoooT
		../systems/FileConverter/converter $4 $trainFile 0 $8 3 $9 ${10}
        fi
	mpiexec $7 ../systems/ParallelRandomForest/exec -trainfile $trainFile -testfile $inputFile -resultfile $5
	
fi
