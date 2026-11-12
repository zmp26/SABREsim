//shows sigmas along y axis
void sigmas_y(TH2* h2){
	gStyle->SetCanvasPreferGL(true);
	//create TBox objects to shade these regions along the bottom axis (or across the plot)
	//get current pad dimensions for Y coordinates of the boxes
	double xmin = h2->GetXaxis()->GetXmin();
	double xmax = h2->GetXaxis()->GetXmax();

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
	//5Li 3/2-
	double sig1_min = 21.7818e6;
	double sig1_max = 21.7915e6;

	double sig2_min = 21.7769e6;
	double sig2_max = 21.7964e6;

	double sig3_min = 21.7721e6;
	double sig3_max = 21.8013e6;

	h2->Draw("COLZ");

	// ±3 sigma background band (Lightest)
	TBox *box3s = new TBox(xmin, sig3_min, xmax, sig3_max);
	box3s->SetFillColorAlpha(color3, 0.2);
	box3s->Draw("SAME");

	// ±2 sigma background band (Medium)
	TBox *box2s = new TBox(xmin, sig2_min, xmax, sig2_max);
	box2s->SetFillColorAlpha(color2, 0.2);
	box2s->Draw("SAME");

	// ±1 sigma background band (Darkest)
	TBox *box1s = new TBox(xmin, sig1_min, xmax, sig1_max);
	box1s->SetFillColorAlpha(color1, 0.2);
	box1s->Draw("SAME");

	//draw the histogram last since i can't figure out how to make the boxes appear transparent and they block the data
	h2->Draw("COLZ SAME");
}