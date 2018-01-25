pipeline {
    agent any

    stages {
	    stage('Building Stage'){
	    	    steps{
		    echo 'Cloning Subrepo Dependencies'
		    dir('LibChemist'){
		    git credentialsId: 'c2b379ef-2e89-4f24-bc4f-2f4d1f85f6d1', url: 'https://github.com/NWChemEx-Project/LibChemist.git', branch: 'master'
		    }
		    
		    sh '''
		    set +x
		    source /etc/profile
		    module load gcc/7.1.0
		    module load cmake

		    echo 'Building LibChemist'
		    mkdir -p root		    
		    cd LibChemist
		    cmake -DCMAKE_INSTALL_PREFIX=../root -H. -Bbuild
		    cd build
		    make && make install
		    cd ../../

		    echo 'Building IntegralsEx'
		    cmake -DCMAKE_PREFIX_PATH=${PWD}/root/usr/local -H. -Bbuild
		    cd build
		    make
		    '''
		    }
		    }

	    stage('Testing Stage'){
		    steps{
		    sh '''
		    set +x
		    source /etc/profile
		    module load cmake
		    
		    cd build
		    ctest
		    '''
		    }
		    }
		    }
		    }
		    
		    

