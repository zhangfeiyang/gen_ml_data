#include "rootheader.h"

double *get_gamma_sigma_mean(TString filename,double a,double b,double p0,double p1,double p2,double p3,double p4,double scale){
//	double scale = 1.30566e3;	
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
			mean = scale*(p0*abs_e + p1*abs_e*abs_e+p2)/(1+p3*exp(-p4*abs_e))*abs_e/E_e;
			//cout << (p0*abs_e + p1*abs_e*abs_e+p2)/(1+p3*exp(-p4*abs_e))*abs_e/E_e << "\t"<<mean<<"\n"; 
			mean = gRandom->PoissonD(mean);
			PE[i] += mean;
		}
		
		double sigma_2 = a*a*PE[i] + b*b*PE[i]*PE[i];
		double sigma = sqrt(sigma_2);
		PE[i] += gRandom->Gaus(PE[i],sigma);
	}
	double rms = TMath::RMS(N,PE);
	double average = TMath::Mean(N,PE);
	double *result = new double[2];
	result[0] = average;
	result[1] = rms;
	f->Close();
	return result;
}
double *get_eplus_sigma_mean(double E,double a,double b,double p0,double p1,double p2,double p3,double p4,double scale){

	double Ek = E - 1.0219978;

	TFile *f = new TFile("/home/zhangfy/gen_ml_data/Ge68.root","read");
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
			mean = scale*(p0*abs_e + p1*abs_e*abs_e+p2)/(1+p3*exp(-p4*abs_e))*abs_e/E_e;
			mean = gRandom->PoissonD(mean);
			double sigma_2 = a*a*mean + b*b*mean*mean;
			double sigma = sqrt(sigma_2);
			PE[i] += gRandom->Gaus(mean,sigma);
		}
		mean = scale*(p0*Ek + p1*Ek*Ek+p2)/(1+p3*exp(-p4*Ek));
		mean = gRandom->PoissonD(mean);
		double sigma_2 = a*a*mean + b*b*mean*mean;
		double sigma = sqrt(sigma_2);
		PE[i] += gRandom->Gaus(mean,sigma);
	}

    double rms = TMath::RMS(N,PE);
    double average = TMath::Mean(N,PE);
    double *result = new double[2];
    result[0] = average;
    result[1] = rms;
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
	TNtuple *t = new TNtuple("t","","a:b:gm0:gm1:gm2:gm3:gm4:gm5:gm6:gm7:egm0:egm1:egm2:egm3:egm4:egm5:egm6:egm7:eplusm2:eplusm3:eplusm4:eplusm5:eplusm6:eplusm7:eplusm8:eeplusm2:eeplusm3:eeplusm4:eeplusm5:eeplusm6:eeplusm7:eeplusm8");
	int N = 1000;
	double mean;
	double emean;
	float *result = new float[32];
	for(int j=0;j<N;j++){
		cout << "finished\t"<<j*100.0/N<<"%\n";	
		double scale = 1200;
		scale = gRandom->Gaus(scale,scale*0.2);
		a = 0.0082509*sqrt(scale);
		b = 0.005691;
		a = gRandom->Gaus(a,a*0.1);		
		b = gRandom->Gaus(b,b*0.1);		
		result[0] = a;
		result[1] = b;
		for(int i=0;i<9;i++){
			data = get_gamma_sigma_mean(filename[i],a,b,par[0],par[1],par[2],par[3],par[4],scale);
			mean = data[0];
			emean = data[1];
			result[i+2] = mean;
			result[i+10] = emean;
		}
		
		for(int i=2;i<9;i++){
            data = get_eplus_sigma_mean(i/1.0,a,b,par[0],par[1],par[2],par[3],par[4],scale);
            mean = data[0];
            emean = data[1];
            result[i+16] = mean;
            result[i+23] = emean;
		}

		t->Fill(result);
	}
	TFile *file = new TFile(Form("/home/zhangfy/gen_ml_data/data2_%i.root",seed),"recreate");
	t->Write();
	file->Close();
	return 0;
}
