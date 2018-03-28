# examples:
# pink: ./cluster.sh ~/homedir/datasets/Credit/credit-card.bin 4 out-pink.txt 4
# kmeans: ./cluster.sh ~/homedir/datasets/Credit/credit-card.bin 4 out-kmeans.txt 3 "-n 10"
# optics: ./cluster.sh ~/homedir/datasets/Credit/credit-card.bin 4 out-optics.txt 2 "-m 5 -e 60"
# dbscan: ./cluster.sh ~/homedir/datasets/Credit/credit-card.bin 4 out.txt 1 "-minpts 2 -e 0.8"

# script parameters are:
# 1- input file
# 2- file type: 0-csv, 1-csv with header, 2-  sparse csv, 3- ARFF file, 4- binary (4-bytes), 5- binary (8-bytes)
# 3- output file
# 4- clustering technique
#    1- dbscan, 2- optics, 3- kmeans, 4- pink, 5-HDBSCAN, 6- Canopy
# 5- clustering parameters

# File types per clustering technique
# dbscan: Binary file; Number of points, N and number of dimensions, D (each 4 bytes). followed by the points coordinates (N x D floating point numbers (4-bytes)).
# optics: Binary file; Number of points, N and number of dimensions, D (each 4 bytes). followed by the points coordinates (N x D floating point numbers (4-bytes)).
# kmeans: Binary file; Number of points, N and number of dimensions, D (each 4 bytes). followed by the points coordinates (N x D floating point numbers (4-bytes)).
# pink:   Binary file; Number of points, N and number of dimensions, D (each 4 bytes). followed by the points coordinates (N x D floating point numbers (8-bytes)).
# HDBSCAN: normal csv file
# Canopy:  normal csv file

inputFile=$1

case $4 in
        1)
            # dbscan
            if [ $2 -ne 4 ]
            then
              echo "File needs conversion ..."
              inputFile=tempooo
              ../systems/FileConverter/converter $1 $inputFile 0 $2 4 -1 -1
            fi
	    executable="../systems/dbscan-v1.0.0/parallel_mpi/mpi_dbscan"
            mpiexec -n 4 $executable -i $inputFile -o $3 -b $5
            ;;
         
        2)
            # optics
            if [ $2 -ne 4 ]
            then
              echo "File needs conversion ..."
              inputFile=tempooo
              ../systems/FileConverter/converter $1 $inputFile 0 $2 4 -1 -1
            fi
	    executable="../systems/optics-v1.0.0/parallel_mpi/mpi_optics"
            mpiexec -n 4 $executable -i $inputFile -o $3 -b $5
            ;;
         
        3)
            # kmeans
            if [ $2 -ne 4 ]
            then
              echo "File needs conversion ..."
              inputFile=tempooo
              ../systems/FileConverter/converter $1 $inputFile 0 $2 4 -1 -1
            fi
            executable="../systems/parallel-kmeans/mpi_main"
            mpiexec -n 4 $executable -i $inputFile -o $3 -b $5
            ;;
        4)
            # pink
            if [ $2 -ne 5 ]
            then
              echo "File needs conversion ..."
              inputFile=tempooo
              ../systems/FileConverter/converter $1 $inputFile 0 $2 5 -1 -1
            fi
            executable="../systems/pink-v1.0.0/pink"
            mpiexec -n 5 $executable $inputFile > temp
            ;;
        5)
            # hdbscan
            if [ $2 -ne 0 ]
            then
              echo "File needs conversion ..."
              inputFile=tempooo
              ../systems/FileConverter/converter $1 $inputFile 0 $2 0 -1 -1
            fi

            executable="../systems/HDBSCAN/hdbscan"
            mpiexec -n 5 $executable -nPartitions 10 -file $inputFile -threads 5 -output $3
            ;;
        6)
            # canopy
            if [ $2 -ne 0 ]
            then
              echo "File needs conversion ..."
              inputFile=tempooo
              ../systems/FileConverter/converter $1 $inputFile 0 $2 0 -1 -1
            fi

            executable="../systems/Canopy/canopy_mpi -file $inputFile $5"
            $executable>$3
            ;;
        *)
            echo "Usage: parameter #3 should be between 1-5"
            exit 1
 
esac


# do post-processing (if needed)
case $4 in
        1)
            # dbscan
            ;;

        2)
            # optics
            ../systems/optics-v1.0.0/mst_to_clusters/mst_to_clusters -i $3 -o temp.txt -t 1 -e 0.5
	    rm $3
            mv temp.txt $3
            ;;

        3)
            # kmeans
	    mv $1.membership $3
            ;;
        4)
            # pink
            ../systems/pink-v1.0.0/thresh k 10 temp $3
            ;;
        5)
            # hdbscan
            ;;
esac
