pipeline {
    agent any

    stages {
	    stage('Checkout'){
		    echo 'Checking out LibChemist'
		    git credentialsId: '422b0eed-700d-444d-961c-1e58cc75cda2'', url: 'https://github.com/NWChemEx-Project/LibChemist.git', branch: 'master'
	    }
  /*      stage('Build') {
            steps {
                echo 'Building...'
	        sh '''
		set +x
	        source /etc/profile
  	        module load gcc/7.1.0-4bgguyp
     	        module load cmake
	        cmake -H. -Bbuild
	        cd build
	        make
	        '''
            }
        }
        stage('Test') {
            steps {
                echo 'Testing..'
	  	sh'''
		set +x
     	        source /etc/profile
	        module load cmake
	        cd build
	        ctest
	        '''*/
            }
        }
    }
}
