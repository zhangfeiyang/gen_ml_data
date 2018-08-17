/// \file
/// \ingroup tutorial_mlp
///  This macro shows the use of an ANN for regression analysis:
/// given a set {i} of input vectors i and a set {o} of output vectors o,
/// one looks for the unknown function f(i)=o.
/// The ANN can approximate this function; TMLPAnalyzer::DrawTruthDeviation
/// methods can be used to evaluate the quality of the approximation.
///
/// For simplicity, we use a known function to create test and training data.
/// In reality this function is usually not known, and the data comes e.g.
/// from measurements.
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Axel Naumann, 2005-02-02


void test3() {
	
	TFile *f = new TFile("data0.root","read");
	TTree *t = (TTree*)f->Get("t");
	
   //TMultiLayerPerceptron* mlp=new TMultiLayerPerceptron("em0/m0,em1/m1,em2/m2,em3/m3,em4/m4,em5/m5,em6/m6,em7/m7:10:8:8:8:a,b",t,"Entry$%2","(Entry$)%2==0");
   TMultiLayerPerceptron* mlp=new TMultiLayerPerceptron("em0/m0,em1/m1,em2/m2,em3/m3,em4/m4,em5/m5,em6/m6,em7/m7,em8/m8:10:8:8:8:8:8:a,b",t,"Entry$%10!=0","(Entry$)%10==0");
   //TMultiLayerPerceptron* mlp=new TMultiLayerPerceptron("x,y,z:10:8:f,g",t,"Entry$%2","(Entry$)%2==0");
   mlp->Train(300,"graph update=10");
	ifstream fin("data0");
	double tmp,sigma,mean,pars[9];
	for(int i=0;i<9;i++){
		fin>>tmp>>mean>>tmp>>sigma>>tmp;
		pars[i] = sigma/mean;
	}
	cout <<mlp->Evaluate(0,pars)<<"\n";
	cout <<mlp->Evaluate(1,pars)<<"\n";
/*	
   double pars[3];
   for (Int_t ix=0; ix<15; ix++) {
	x = -1 + 2./30*(2*ix+1);
	y = -1 + 2./30*(2*ix+1);
	z = -1 + 2./30*(2*ix+1);
	pars[0] = x;
	pars[1] = y;
	pars[2] = z;
	data = theUnknownFunction(x,y,z);
	cout << x <<"\t"<<theUnknownFunction(x,y,z)<<"\t" <<mlp->Evaluate(0,pars)<<"\n";
   }
*/
}
