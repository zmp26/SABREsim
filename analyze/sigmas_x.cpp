//shows sigmas along x axis
void sigmas_x(TH2* h2){
	gStyle->SetCanvasPreferGL(true);
	//create TBox objects to shade these regions along the bottom axis (or across the plot)
	//get current pad dimensions for Y coordinates of the boxes
	double ymin = h2->GetYaxis()->GetXmin();
	double ymax = h2->GetYaxis()->GetXmax();

	//green:
	// Color_t color1 = TColor::GetColor("#0A2F1D"); // Dark green
	// Color_t color2 = TColor::GetColor("#1B5E20"); // Medium green
	// Color_t color3 = TColor::GetColor("#388E3C"); // Light tint

	//lsu purple:
	//Color_t color3 = TColor::GetColor("#C2AFE3"); // Light purple (Thistle) for ±3 sigma
	//Color_t color2 = TColor::GetColor("#7A52A5"); // Medium purple (Medium Purple) for ±2 sigma
	//Color_t color1 = TColor::GetColor("#461D7C"); // Dark purple (LSU purple) for ±1 sigma

	//lsu gold:
	Color_t color1 = TColor::GetColor("#C69214"); // Dark green
	Color_t color2 = TColor::GetColor("#FDD023"); // Medium green
	Color_t color3 = TColor::GetColor("#FFE680"); // Light tint

	//define upper/lower bounds for 1/2/3 sigma ranges:
	
	//8Be 0+
	/*
	double sig1_min = 55.5746e6;
	double sig1_max = 55.5750e6;

	double sig2_min = 55.5744e6;
	double sig2_max = 55.5752e6;

	double sig3_min = 55.5742e6;
	double sig3_max = 55.5754e6;
	*/

	//8Be 0+
	double sig1_min = 55.610e6;
	double sig1_max = 55.630e6;

	double sig2_min = 55.600e6;
	double sig2_max = 55.639e6;

	double sig3_min = 55.591e6;
	double sig3_max = 55.649e6;
	

	//5Li 3/2-
	/*
	double sig1_min = e6;
	double sig1_max = e6;

	double sig2_min = e6;
	double sig2_max = e6;

	double sig3_min = e6;
	double sig3_max = e6;
	*/

	h2->Draw("COLZ");

	// ±3 sigma background band (Lightest)
	TBox *box3s = new TBox(sig3_min, ymin, sig3_max, ymax);
	box3s->SetFillColorAlpha(color3, 0.2);
	box3s->Draw("SAME");

	// ±2 sigma background band (Medium)
	TBox *box2s = new TBox(sig2_min, ymin, sig2_max, ymax);
	box2s->SetFillColorAlpha(color2, 0.2);
	box2s->Draw("SAME");

	// ±1 sigma background band (Darkest)
	TBox *box1s = new TBox(sig1_min, ymin, sig1_max, ymax);
	box1s->SetFillColorAlpha(color1, 0.2);
	box1s->Draw("SAME");

	//draw the histogram last since i can't figure out how to make the boxes appear transparent and they block the data
	h2->Draw("COLZ SAME");
}