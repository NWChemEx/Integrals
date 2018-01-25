pipeline {
    agent any

    stages {
	    stage('Checkout'){
	    	    steps{
		    echo 'Cloning and Building Subrepo Dependencies'
		    git credentialsId: '422b0eed-700d-444d-961c-1e58cc75cda2', url: 'https://github.com/NWChemEx-Project/LibChemist.git', branch: 'master'
		    sh '''
		    set +x
		    source /etc/profile
		    module load gcc/7.1.0
		    module load cmake
		    mkdir -p root
		    cd LibChemist
		    cmake -DCMAKE_INSTALL_PREFIX=../root -H. -Bbuild
		    cd build
		    make && make install
		    '''
		    }
		    }
	    stage('Build'){
		    steps{
		    echo 'Building...'
		    sh '''
		    echo $pwd
		    ls
		    '''
		    }
		    }
		    }
		    }
		    
		    

