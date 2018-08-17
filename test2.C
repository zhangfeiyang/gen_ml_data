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

Double_t theUnknownFunction(Double_t x,double y,double z) {
   return (sin(x)*y+z);
}

void test2() {
   // create a tree with train and test data.
   // we have two input parameters x and y,
   // and one output value f(x,y)
   TNtuple* t=new TNtuple("tree","tree","x:y:z:f");
   TRandom r;
   for (Int_t i=0; i<10000; i++) {
      Float_t x=r.Uniform(-1,1);
      Float_t y=r.Uniform(-1,1);
      Float_t z=r.Uniform(-1,1);
      t->Fill(x,y,z,theUnknownFunction(x,y,z));
   }

   TMultiLayerPerceptron* mlp=new TMultiLayerPerceptron("x,y,z:10:8:f",t,
      "Entry$%2","(Entry$)%2==0");
   mlp->Train(150,"graph update=10");
	
	double pars[3];
   for (Int_t ix=0; ix<15; ix++) {
	double x = -1 + 2./30*(2*ix+1);
	double y = -1 + 2./30*(2*ix+1);
	double z = -1 + 2./30*(2*ix+1);
	pars[0] = x;
	pars[1] = y;
	pars[2] = z;
	cout << x <<"\t"<<theUnknownFunction(x,y,z)<<"\t" <<mlp->Evaluate(0,pars)<<"\n";
   }
}
