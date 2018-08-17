#include "rootheader.h"

double *get_sigma_mean(TString filename,double a,double b,double p0,double p1,double p2,double p3,double p4){
	double scale = 1.30566e3;	
	TFile *f = new TFile(filename,"read");
	TTree *t = (TTree*)f->Get("Gamma2e");
	double mean;
	int N = 10000;
	double PE[N];
	for(int i=0;i<N;i++){
		t->GetEntry(i);
		PE[i] = 0;
		int nEle = t->GetLeaf("nEle")->GetValue(0);
		for(int j=0;j<nEle;j++){
			double E_e = t->GetLeaf("E_e")->GetValue(j);
			double abs_e = TMath::Abs(E_e);
			if(abs_e==0) continue;
			mean = (p0*abs_e + p1*abs_e*abs_e+p2)/(1+p3*exp(-p4*abs_e))*abs_e/E_e;
			abs_e = mean;
			double E0 = 1259.53/1305.66;
			double sigma_2 = a*a*(abs_e+E0)+b*b*(abs_e+E0)*(abs_e+E0)-37.6426*37.6426;
			double sigma = sqrt(sigma_2);

			PE[i] += gRandom->Gaus(mean,sigma);
			PE[i] += gRandom->Gaus(1259.53,37.6426);
		}
	}
	double rms = TMath::RMS(N,PE);
	double average = TMath::Mean(N,PE);
	double *result = new double[2];
	result[0] = average*scale;
	result[1] = rms*scale;
	f->Close();
	return result;
}

//void gen_data(){
int main(int argc,char **argv){

	int seed = atoi(argv[1]);
	gRandom->SetSeed(seed);

	ifstream finfile("/home/zhangfy/gen_ml_data/log");
	TString filename[9];
	for(int i=0;i<9;i++) finfile>>filename[i];

	//double par[5] = {1.06269,0.00444552,-0.00150689,0.0965178,1.8072};
	double par[5] = {1.071,2.393e-3,-1.765e-3,8.816e-2,1.402};
	double a = 0.0282509;
	double b = 0.005691;
	double *data;
	TNtuple *t = new TNtuple("t","","a:b:m0:m1:m2:m3:m4:m5:m6:m7:m8:em0:em1:em2:em3:em4:em5:em6:em7:em8");
	int N = 1000;
	double mean[9];
	double emean[9];
	float *result = new float[20];
	for(int j=0;j<N;j++){
		cout << "finished\t"<<j*100.0/N<<"%\n";	
		a = 0.0282509;
		b = 0.005691;
		a = gRandom->Gaus(a,a*0.1);		
		b = gRandom->Gaus(b,b*0.1);		
		result[0] = a;
		result[1] = b;
		for(int i=0;i<9;i++){
			cout << filename[i] <<"\n";
			data = get_sigma_mean(filename[i],a,b,par[0],par[1],par[2],par[3],par[4]);
			mean[i] = data[0];
			emean[i] = data[1];
			result[i+2] = mean[i];
			result[i+11] = emean[i];
		//	cout<<data[0]<<"\t"<<data[1]<<"\n";
		}
		t->Fill(result);
	}
	TFile *file = new TFile(Form("/home/zhangfy/gen_ml_data/data1_%i.root",seed),"recreate");
	t->Write();
	file->Close();
	return 0;
}
