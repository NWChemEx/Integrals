pipeline {
    agent any

    stages {
	    stage('Checkout'){
		    echo 'Checking out LibChemist'
		    git credentialsId: '422b0eed-700d-444d-961c-1e58cc75cda2'', url: 'https://github.com/NWChemEx-Project/LibChemist.git', branch: 'master'
			    }
    }
}
