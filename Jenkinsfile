pipeline {
    agent any

    stages {
	    stage('Building Stage'){
	    	    steps{
		    echo 'Cloning Subrepo Dependencies'
		    dir('LibChemist'){
		    git credentialsId: '422b0eed-700d-444d-961c-1e58cc75cda2', url: 'https://github.com/NWChemEx-Project/LibChemist.git', branch: 'master'
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
		    
		    

