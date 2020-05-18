#include<math.h>
#include<vector>
#include<iostream>
#include<string>
#include<fstream>
#include<chrono>


using namespace std;


void PrintBanner(){
    printf("\n"
"                                                                \n"
"           ____ ___.____       _________________                \n" 
"          |    |   \\    |     /   _____/\\______ \\            \n" 
"  ______  |    |   /    |     \\_____  \\  |    |  \\    ______ \n" 
" /_____/  |    |  /|    |___  /        \\ |    `   \\  /_____/  \n" 
"          |______/ |_______ \\/_______  //_______  /            \n" 
"                           \\/        \\/         \\/           \n" 
"                                                                \n"  
"\n  Unfolding Level Spacing Distributin Calculator\n\n");
}


const double S_CRIT=0.508;

struct arg_list{
    int num_eigval = 100;
    int spectral_bins = 30;
    int num_resampling = 50;
    string infile = "";
    string outfile = "";    

    friend ostream& operator<<(ostream& os, const arg_list& al);
};


ostream& operator<<(ostream& os, const arg_list& al){
    os<<"Number of eigenvalues: "<<al.num_eigval<<endl;
    os<<"Number of spectral bins: "<<al.spectral_bins<<endl;
    os<<"Number of resamplings: "<<al.num_resampling<<endl;
    os<<"Input datafile: "<<al.infile<<endl;
    os<<"Output datafile: "<<al.outfile<<endl;
    return os;
}

void parse_arguments(arg_list& arg, char** argv){
    //int fixed_args = 4;
    arg.num_eigval = atoi(argv[1]);
    arg.spectral_bins = atoi(argv[2]);
    arg.num_resampling = atoi(argv[3]);
    arg.infile = argv[4];
    arg.outfile = argv[5];

    //Arguments checkings:
    if(arg.infile == ""){
        throw "ERROR: inputfile is empty!";
    }

    if(arg.outfile == ""){
        throw "ERROR: outputfile is empty!";
    }

    if(arg.num_eigval <= 0){
        throw "ERROR: Negative number of eigenvalues!";
    }

    if(arg.spectral_bins <= 0){
        throw "ERROR: Negative number of spectral bins!";
    }

}

// This function Load the input data
void LoadInput(arg_list arg, vector<double>& total_eigs ){
    FILE* input_eigs = fopen(arg.infile.c_str(), "r");
    double aux_eig;

    if(input_eigs == NULL){
        printf("ERROR! Could not open the input file %s\n\n", arg.infile.c_str());
        exit(1);
    }    

    while((fscanf(input_eigs, "%lf", &aux_eig) == 1)){
        total_eigs.push_back(aux_eig);
    }    

    fclose(input_eigs);

    if(total_eigs.size()%arg.num_eigval != 0){
        printf("ERROR: Total eigs can not be divided by num eigs!");
        exit(1);
    }
}


// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(vector<double>& arr, vector<int>& conf_array, int l, int m, int r){ 
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 =  r - m; 
  
    /* create temp arrays */
    double L[n1], R[n2], CL[n1], CR[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){ 
        L[i] = arr[l + i]; 
        CL[i] = conf_array[l + i]; 
    }
    for (j = 0; j < n2; j++){ 
        R[j] = arr[m + 1+ j];
        CR[j] = conf_array[m + 1+ j];
    } 
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2){ 
        if (L[i] <= R[j]){ 
            arr[k] = L[i]; 
            conf_array[k] = CL[i]; 
            i++; 
        } 
        else{ 
            arr[k] = R[j]; 
            conf_array[k] = CR[j]; 
            j++; 
        } 
        k++; 
    } 
  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1){ 
        arr[k] = L[i]; 
        conf_array[k] = CL[i]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2){ 
        arr[k] = R[j]; 
        conf_array[k] = CR[j]; 
        j++; 
        k++; 
    } 
} 



/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort(vector<double>& arr, vector<int>& conf_array, int l, int r){ 
    if (l < r){ 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        int m = l+(r-l)/2; 
  
        // Sort first and second halves 
        mergeSort(arr, conf_array, l, m); 
        mergeSort(arr, conf_array, m+1, r); 
  
        merge(arr, conf_array, l, m, r); 
    } 
} 

// This function initialize the array with the configuration index
void InitConfArray(vector<int>& conf_array, arg_list arg, int len){
    for(int i=0;i<len; i++){
        conf_array[i]=(i/arg.num_eigval);
    }
}

//This function creates an array with the spectral bins using total_eigsi
//total_eigs array must be ordered!!
//TODO: add a check to see if total_eigs is ordered or not
void InitBinArray(vector<double>& bins, vector<double>& total_eigs, arg_list arg){
    double min_eig = total_eigs[0];
    double max_eig = total_eigs[total_eigs.size()-1];
    double size_block = (max_eig - min_eig)/arg.spectral_bins;

    for(unsigned int i=0;i<bins.size();i++){
        bins[i] = min_eig + size_block*i;
    }
}

//This function saves the ranks in the matrix ulsd using the fixed bin [sx_bin:dx_bin]
//TODO: ADD SOME CHECKS!
void ComputeUlsdFixedBin(vector<vector<double>>& conf_rank, vector<double>& total_eigs, 
                        vector<int>& conf_array, int num_conf, double sx_bin, 
                        double dx_bin, vector<double>& ulsd, int* spec_counter){
    
    *spec_counter = 0;
    //find the first eigenvalue in the bin
    int offset =0;
    while(total_eigs[offset]<sx_bin){
        offset+=1;
    }

    //save the ranks for each conf
    int index=offset;
    while(total_eigs[index]>sx_bin && total_eigs[index]<dx_bin){
        conf_rank[conf_array[index]].push_back((double)(index+1)/num_conf);
        index+=1;
        *spec_counter += 1;
    }
 
	//add the extremal lambda
	vector<int> stopper;
	int aux;
	int condition;	

	while((stopper.size()< (unsigned int) num_conf) && ((unsigned int)index<total_eigs.size()-1)){
		aux = conf_array[index];
		condition=0;
		for(unsigned int k=0; k<stopper.size();k++){
			if(aux==stopper[k]){
				condition=1;
			}	
		}
		
		if(condition==1){
			index+=1;
		}
		else{
			conf_rank[conf_array[index]].push_back((double)(index+1)/num_conf);
			stopper.push_back(aux);
			index +=1;
		}

	}
		

    //compute the ulsd
    for(int i=0;i<num_conf;i++){
        if(conf_rank[i].size()>1){
            for(unsigned int j=0;j<(conf_rank[i].size()-1);j++){
                ulsd.push_back(conf_rank[i][j+1]-conf_rank[i][j]);
            }
        }
    }
}

// This function prints the matrix conf_array, the name ulsd is measleading (for a fixed spectral bin)
void PrintUlsdFixedBin(vector<vector<double>>& ulsd, int num_conf){
    for(int j=0;j<num_conf;j++){
        printf("Conf %d: ", j);
        for(unsigned int i=0;i<ulsd[j].size();i++){
           printf("%lf ", ulsd[j][i]);
        }
       printf("\n");
    } 
}


double Is0Calculator(vector<double>& ulsd){
    double counter=0.0;
    for(unsigned int i=0;i<ulsd.size();i++){
        if(ulsd[i]<S_CRIT){
            counter+=(double)1/ulsd.size();
        }
    }
    return counter;
}


void CreateSampleArray(vector<double>&total_eigs, vector<double>& sampled_eigs, arg_list arg, int num_conf){
    int random;
    for(int i=0;i<num_conf;i++){
        random = rand()%num_conf;
        for(int j=0; j<arg.num_eigval;j++){
            sampled_eigs[i*arg.num_eigval + j] = total_eigs[random*arg.num_eigval +j];
        }   
    }
}


void Is0Boostrap(vector<double>& total_eigs, arg_list args, int num_conf, 
                 int num_resampling, 
                 vector<double>& Is0, vector<double>& err_Is0, 
                 vector<double>& spec_den, vector<double>& err_spec_den, 
 				 vector<double>& ave_s, vector<double>& err_ave_s,
                 vector<double>& array_bin){

    //Initialize stuff
    vector<double> sampled_eigs(total_eigs.size(), 0.0);
    vector<double> bins(args.spectral_bins+1, 0.0);
    vector<double> aux(total_eigs.size());
    vector<double> ulsd;
    vector<double> Is0_aux(bins.size()-1, 0.0);
    vector<double> errIs0_aux(bins.size()-1, 0.0); 
    vector<double> spec_den_aux(bins.size()-1, 0.0);
    vector<double> err_spec_den_aux(bins.size()-1, 0.0); 
    vector<double> ave_s_aux(bins.size()-1, 0.0);
    vector<double> err_ave_s_aux(bins.size()-1, 0.0); 
    vector<int> conf_array(total_eigs.size(), 0);
    vector<vector<double>> conf_rank(num_conf);
    int spec_den_counter;
    double aux2;
    double interval;

    //Create the array with the bins
    for(unsigned int a =0; a<total_eigs.size();a++){
        aux[a] = total_eigs[a];
    }
    mergeSort(aux, conf_array, 0, total_eigs.size()-1);
    InitBinArray(bins, aux, args);     
       
    for(unsigned int a=0;a<bins.size()-1;a++){
        array_bin.push_back((bins[a+1]+bins[a])/2);
    }

    interval = bins[2] - bins[1];
    
    for(int res=0;res<num_resampling;res++){
        CreateSampleArray(total_eigs, sampled_eigs, args, num_conf);

        //Initialize the confs array
        InitConfArray(conf_array, args, total_eigs.size());
        
        //Sort the eigenvalues
        mergeSort(sampled_eigs, conf_array, 0, total_eigs.size()-1);

        //cycle the bins
        for(unsigned int bbin=0;bbin<bins.size()-1;bbin++){
            ComputeUlsdFixedBin(conf_rank, sampled_eigs, conf_array, num_conf, bins[bbin], bins[bbin+1], ulsd, &spec_den_counter);
           
            //Divide by the mean
            aux2=0.0;
            for(unsigned int b =0;b<ulsd.size();b++){
               aux2 = aux2 + ulsd[b]/ulsd.size();
            }
           
 
            //actual computation of Is0
			ave_s_aux[bbin] += aux2/num_resampling;
			err_ave_s_aux[bbin] += pow(aux2, 2 );
            Is0_aux[bbin]+=((Is0Calculator(ulsd))/num_resampling);
            errIs0_aux[bbin]+=pow((Is0Calculator(ulsd)), 2);
            spec_den_aux[bbin] += (double) spec_den_counter/num_resampling; 
            err_spec_den_aux[bbin] += pow( (double) spec_den_counter, 2); 

            for(int i=0; i<num_conf;i++){
                conf_rank[i].resize(0);
            }
            ulsd.resize(0);
        }
    }

    for(unsigned int i=0;i<Is0_aux.size();i++){
     	ave_s.push_back(ave_s_aux[i]);
		err_ave_s.push_back(sqrt(err_ave_s_aux[i]/num_resampling - pow(ave_s_aux[i],2)));
	    Is0.push_back(Is0_aux[i]);
        err_Is0.push_back(sqrt(errIs0_aux[i]/num_resampling - pow(Is0_aux[i],2)));
        spec_den.push_back(spec_den_aux[i]/interval);
        err_spec_den.push_back(sqrt(err_spec_den_aux[i]/num_resampling - pow(spec_den_aux[i],2))/interval);
    }


}


arg_list args;


int main(int argc, char** argv){

    if(argc<6){
        printf("Usage: %s <num eigvals> <num spectral bins> <num num_resampling> <inputfile> <outputfile>\n\n", argv[0]);
        exit(1);
    }
    
    //Initialize the random number generator
    srand (time(NULL));

    //parse the arguments
    parse_arguments(args, argv);

    int num_eigval = args.num_eigval;
    string  otufilename(args.outfile);
    string  infile(args.infile);

    //Variables
    vector<double> total_eigs;
    vector<double> Is0;
    vector<double> err_Is0;
    vector<double> spec_den;
    vector<double> err_spec_den;
	vector<double> ave_s;
	vector<double>err_ave_s;
    vector<double> array_bins;
    int num_conf;
    
    //Banner and arguments
    PrintBanner();
    cout<<"Arguments:\n"<<args<<endl;

    //take starting time
    auto start = std::chrono::high_resolution_clock::now();

    //Load the eigenvalues
    LoadInput(args, total_eigs);        
    num_conf = total_eigs.size()/num_eigval;
    cout<<"Total number of conf: "<<total_eigs.size()/args.num_eigval<<endl;

    Is0Boostrap(total_eigs, args, num_conf, args.num_resampling, Is0, err_Is0, spec_den, err_spec_den, ave_s, err_ave_s, array_bins);

    FILE* out = fopen(args.outfile.c_str(), "w");

    if(out == NULL){
        printf("ERROR! Could not open the output file %s\n\n", args.outfile.c_str());
        exit(1);
    }    


    for(unsigned int i=0;i<Is0.size();i++){
        fprintf(out, "%.16lg %.16lg %.16lg %.16lg %.16lg %.16lg %.16lg\n", array_bins[i], Is0[i], err_Is0[i], spec_den[i]/num_conf, err_spec_den[i]/num_conf, ave_s[i], err_ave_s[i]);
    }

    fclose(out);

     auto end = std::chrono::high_resolution_clock::now();

     std::chrono::duration<double> diff = end-start;

     cout<<"Total time: "<<diff.count()<< "sec"<<endl;

    cout<<"     All fine!";
    cout<<"     ULSD!"<<endl;





}


