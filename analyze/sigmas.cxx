void sigmas(TH2* h2){
	gStyle->SetCanvasPreferGL(true);
	//create TBox objects to shade these regions along the bottom axis (or across the plot)
	//get current pad dimensions for Y coordinates of the boxes
	double ymin = h2->GetYaxis()->GetXmin();
	double ymax = h2->GetYaxis()->GetXmax();
	double xmin = h2->GetXaxis()->GetXmin();
	double xmax = h2->GetXaxis()->GetXmax();

	//green:
	Color_t green1 = TColor::GetColor("#0A2F1D"); // Dark green
	Color_t green2 = TColor::GetColor("#1B5E20"); // Medium green
	Color_t green3 = TColor::GetColor("#388E3C"); // Light tint

	//lsu purple:
	Color_t purple3 = TColor::GetColor("#C2AFE3"); // Light purple (Thistle) for ±3 sigma
	Color_t purple2 = TColor::GetColor("#7A52A5"); // Medium purple (Medium Purple) for ±2 sigma
	Color_t purple1 = TColor::GetColor("#461D7C"); // Dark purple (LSU purple) for ±1 sigma

	//lsu gold:
	Color_t gold1 = TColor::GetColor("#C69214"); // Dark green
	Color_t gold2 = TColor::GetColor("#FDD023"); // Medium green
	Color_t gold3 = TColor::GetColor("#FFE680"); // Light tint

	//define upper/lower bounds for 1/2/3 sigma ranges:
	
	//8Be 0+
	double sig1_8Begs_min = 55.5746e6;
	double sig1_8Begs_max = 55.5750e6;

	double sig2_8Begs_min = 55.5744e6;
	double sig2_8Begs_max = 55.5752e6;

	double sig3_8Begs_min = 55.5742e6;
	double sig3_8Begs_max = 55.5754e6;
	

	//8Be 2+
	double sig1_8Be2plus_min = 55.610e6;
	double sig1_8Be2plus_max = 55.630e6;

	double sig2_8Be2plus_min = 55.600e6;
	double sig2_8Be2plus_max = 55.639e6;

	double sig3_8Be2plus_min = 55.591e6;
	double sig3_8Be2plus_max = 55.649e6;
	

	//5Li 3/2-
	double sig1_5Ligs_min = 21.7818e6;
	double sig1_5Ligs_max = 21.7915e6;

	double sig2_5Ligs_min = 21.7769e6;
	double sig2_5Ligs_max = 21.7964e6;

	double sig3_5Ligs_min = 21.7721e6;
	double sig3_5Ligs_max = 21.8013e6;

	h2->Draw("COLZ");

	//8Begs:
	// ±3 sigma background band (Lightest)
	TBox *box3s_8Begs = new TBox(sig3_8Begs_min, ymin, sig3_8Begs_max, ymax);
	box3s_8Begs->SetFillColorAlpha(purple3, 0.2);
	box3s_8Begs->Draw("SAME");

	// ±2 sigma background band (Medium)
	TBox *box2s_8Begs = new TBox(sig2_8Begs_min, ymin, sig2_8Begs_max, ymax);
	box2s_8Begs->SetFillColorAlpha(purple2, 0.2);
	box2s_8Begs->Draw("SAME");

	// ±1 sigma background band (Darkest)
	TBox *box1s_8Begs = new TBox(sig1_8Begs_min, ymin, sig1_8Begs_max, ymax);
	box1s_8Begs->SetFillColorAlpha(purple1, 0.2);
	box1s_8Begs->Draw("SAME");

	//8Be 2+:
	// ±3 sigma background band (Lightest)
	TBox *box3s_8Be2plus = new TBox(sig3_8Be2plus_min, ymin, sig3_8Be2plus_max, ymax);
	box3s_8Be2plus->SetFillColorAlpha(green3, 0.2);
	box3s_8Be2plus->Draw("SAME");

	// ±2 sigma background band (Medium)
	TBox *box2s_8Be2plus = new TBox(sig2_8Be2plus_min, ymin, sig2_8Be2plus_max, ymax);
	box2s_8Be2plus->SetFillColorAlpha(green2, 0.2);
	box2s_8Be2plus->Draw("SAME");

	// ±1 sigma background band (Darkest)
	TBox *box1s_8Be2plus = new TBox(sig1_8Be2plus_min, ymin, sig1_8Be2plus_max, ymax);
	box1s_8Be2plus->SetFillColorAlpha(green1, 0.2);
	box1s_8Be2plus->Draw("SAME");

	//5Ligs:
	// ±3 sigma background band (Lightest)
	TBox *box3s_5Ligs = new TBox(xmin, sig3_5Ligs_min, xmax, sig3_5Ligs_max);
	box3s_5Ligs->SetFillColorAlpha(gold3, 0.2);
	box3s_5Ligs->Draw("SAME");

	// ±2 sigma background band (Medium)
	TBox *box2s_5Ligs = new TBox(xmin, sig2_5Ligs_min, xmax, sig2_5Ligs_max);
	box2s_5Ligs->SetFillColorAlpha(gold2, 0.2);
	box2s_5Ligs->Draw("SAME");

	// ±1 sigma background band (Darkest)
	TBox *box1s_5Ligs = new TBox(xmin, sig1_5Ligs_min, xmax, sig1_5Ligs_max);
	box1s_5Ligs->SetFillColorAlpha(gold1, 0.2);
	box1s_5Ligs->Draw("SAME");

	//draw the histogram last since i can't figure out how to make the boxes appear transparent and they block the data
	//h2->Draw("COLZ SAME");
}