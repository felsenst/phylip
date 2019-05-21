package phylip;

import java.awt.EventQueue;
import javax.swing.JFileChooser;
import java.io.File;
import javax.swing.JFrame;
import javax.swing.JButton;
import javax.swing.JTextField;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JOptionPane;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import javax.swing.JComboBox;
import javax.swing.SwingUtilities;

import java.awt.Font;
import javax.swing.SwingConstants;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import com.sun.jna.Library;
import com.sun.jna.Native;

import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class DnaMLKUserInterface {
   public interface DnaMLK extends Library {
        public void dnamlk(
        		String infile,
        		String intree,
        		String wgtsfile,
        		String catsfile,
        		String outfile,
        		String outfileopt,
        		String outtree,
        		String outtreeopt,
        		String TreeUseMethod,
        		boolean UseLengths,
        		double TTratio,
        		boolean useEmpBF,
        		double BaseFreqA,
        		double BaseFreqC,
        		double BaseFreqG,
        		double BaseFreqTU,
        		boolean OneCat,
        		int NumCats,
        		// these are explicitly named because JNA doesn't pass arrays gracefully
        		double SiteRate1,  
        		double SiteRate2,  
        		double SiteRate3,  
        		double SiteRate4,  
        		double SiteRate5,  
        		double SiteRate6,  
        		double SiteRate7,  
        		double SiteRate8,  
        		double SiteRate9,  
        		//
        		String RateVar,
        		boolean AdjCor,
        		double BlockLen,
        		double CoeffVar,
        		int NumRates,
        		// these are explicitly named because JNA doesn't pass arrays gracefully
        		double HMMrate1,
                double HMMrate2,
                double HMMrate3,
                double HMMrate4,
                double HMMrate5,
                double HMMrate6,
                double HMMrate7,
                double HMMrate8,
                double HMMrate9,
                //
        		// these are explicitly named because JNA doesn't pass arrays gracefully
                double HMMprob1,
                double HMMprob2,
                double HMMprob3,
                double HMMprob4,
                double HMMprob5,
                double HMMprob6,
                double HMMprob7,
                double HMMprob8,
                double HMMprob9,
                //
                double InvarFract,
        		boolean SitesWeight,
        		boolean GlobalRe,
        		boolean RandInput,
        		int RandNum,
        		int Njumble,
        		boolean MultData,
        		boolean MultDSet,
        		int NumSeqs,
        		boolean InputSeq,
        		boolean PrintData,
        		boolean PrintInd,
        		boolean PrintTree,
        		boolean DotDiff,
        		boolean WriteTree,
        		boolean RecHypo);
    }

	public class DnaMLKData {
		String infile;
		String intree;
		String wgtsfile;
		String catsfile;
		String outfile;
		String outfileopt;
		String outtree;
		String outtreeopt;
		String TreeUseMethod;
		boolean UseLengths;
		double TTratio;
		boolean useEmpBF;
		double BaseFreqA;
		double BaseFreqC;
		double BaseFreqG;
		double BaseFreqTU;
		boolean OneCat;
		int NumCats;
		
		// these are explicitly named because JNA doesn't pass arrays gracefully
		double SiteRate1;  
		double SiteRate2;  
		double SiteRate3;  
		double SiteRate4;  
		double SiteRate5;  
		double SiteRate6;  
		double SiteRate7;  
		double SiteRate8;  
		double SiteRate9;  
		//
		
		String RateVar;
		boolean AdjCor;
		double BlockLen;
		double CoeffVar;
		int NumRates;
		
		// these are explicitly named because JNA doesn't pass arrays gracefully
		double HMMrate1;
        double HMMrate2;
        double HMMrate3;
        double HMMrate4;
        double HMMrate5;
        double HMMrate6;
        double HMMrate7;
        double HMMrate8;
        double HMMrate9;
		//
        
		// these are explicitly named because JNA doesn't pass arrays gracefully
		double HMMprob1;
        double HMMprob2;
        double HMMprob3;
        double HMMprob4;
        double HMMprob5;
        double HMMprob6;
        double HMMprob7;
        double HMMprob8;
        double HMMprob9;
        //
        
        double InvarFract;
		boolean SitesWeight;
		boolean GlobalRe;
		boolean RandInput;
		int RandNum;
		int Njumble;
		boolean MultData;
		boolean MultDSet;
		int NumSeqs;
		boolean InputSeq;
		boolean PrintData;
		boolean PrintInd;
		boolean PrintTree;
		boolean DotDiff;
		boolean WriteTree;
		boolean RecHypo;
	}

	private DnaMLKData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean ExplicitWgts;
	private boolean phylipCall;

	private double[] SiteRateValues;
	private double[] HMMRateValues;
	private double[] HMMProbValues;
	private int lastsiteratecat;
	private int lasthmmcat;

	private JFrame frmDnaMLKControls;
	private JButton btnInputFile;
	private JTextField txtInputFile;
	private JButton btnOutputFile;
	private JTextField txtOutputFile;
	private JLabel lblSearchBest;
	private JLabel lblTTratio;
	private JTextField txtTTratio;
	private JLabel lblEmpBF;
	private JRadioButton rdbtnEmpBFYes;
	private JRadioButton rdbtnEmpBFNo;
	private JLabel lblBaseFreq;
	private JLabel lblBaseFreqA;
	private JTextField txtBaseFreqA;
	private JLabel lblBaseFreqC;
	private JTextField txtBaseFreqC;
	private JLabel lblBaseFreqG;
	private JTextField txtBaseFreqG;
	private JLabel lblBaseFreqTU;
	private JTextField txtBaseFreqTU;
	private JLabel lblCatSites;
	private JRadioButton rdbtnCatYes;
	private JRadioButton rdbtnCatNo;
	private JLabel lblNumCat;
	private JComboBox cmbxNumCat;
	private JLabel lblRateSite;
	private JComboBox cmbxRateSite;
	private JLabel lblSitesWeight;
	private JRadioButton rdbtnSitesWeightYes;
	private JRadioButton rdbtnSitesWeightNo;
	private JLabel lblGlobalRe;
	private JRadioButton rdbtnGlobalReYes;
	private JRadioButton rdbtnGlobalReNo;
	private JLabel lblRandOrder;
	private JRadioButton rdbtnRandOrderYes;
	private JRadioButton rdbtnRandOrderNo;
	private JLabel lblRandSeed;
	private JTextField txtRandSeed;
	private JLabel lblAnalyzeMultData;
	private JRadioButton rdbtnAnalyzeMultDataYes;
	private JRadioButton rdbtnAnalyzeMultDataNo;
	private JLabel lblMultData;
	private JRadioButton rdbtnDataSets;
	private JRadioButton rdbtnWeights;
	private JLabel lblInputSeq;
	private JRadioButton rdbtnInputSeqYes;
	private JRadioButton rdbtnInputSeqNo;
	private JLabel lblPrintData;
	private JRadioButton rdbtnPrintDataYes;
	private JRadioButton rdbtnPrintDataNo;
	private JLabel lblPrintInd;
	private JRadioButton rdbtnPrintIndYes;
	private JRadioButton rdbtnPrintIndNo;
	private JLabel lblPrintTree;
	private JRadioButton rdbtnPrintTreeYes;
	private JRadioButton rdbtnPrintTreeNo;
	private JLabel lblWriteTree;
	private JRadioButton rdbtnWriteTreeYes;
	private JRadioButton rdbtnWriteTreeNo;
	private JLabel lblRecHypo;
	private JRadioButton rdbtnRecHypoYes;
	private JRadioButton rdbtnRecHypoNo;
	private JButton btnExecute;
	private JButton btnQuit;
	private JButton btnWeightFile;
	private JTextField txtWeightFile;
	private JTextField txtInputTree;
	private JButton btnInputTree;
	private JButton btnOutputTree;
	private JTextField txtOutTree;
	private JButton btnCatFile;
	private JTextField txtCatFile;
	private JLabel lblRate;
	private JTextField txtSiteRate;
	private JLabel lblcatnum;
	private JComboBox ratecatnum;
	private JComboBox TreeSearchMethod;
	private JLabel lblHowManyData;
	private JTextField txtNumSeqs;
	private JLabel lblNumberJumble;
	private JTextField txtNumberJumble;
	private JLabel lblRatesAdjCor;
	private JRadioButton rdbtnAdjCorYes;
	private JRadioButton rdbtnAdjCorNo;
	private JLabel lblMeanLen;
	private JTextField txtBlockLen;
	private JLabel lblDotDiff;
	private JRadioButton rdbtnDotDiffYes;
	private JRadioButton rdbtnDotDiffNo;
	private JLabel lblUseLengths;
	private JRadioButton rdbtnUseLengthsYes;
	private JRadioButton rdbtnUseLengthsNo;
	private JLabel lblCoeffVar;
	private JTextField txtCoeffVar;
	private JLabel lblCoeffVarNote;
	private JLabel lblHMM;
	private JLabel lblInvar;
	private JLabel lblNumHMMcats;
	private JComboBox cmbxHMMcount;
	private JLabel lblHMMCat;
	private JComboBox cmbxHMMCat;
	private JLabel lblHMMRate;
	private JTextField txtHMMRate;
	private JLabel lblFracInvar;
	private JTextField txtFracInvar;
	private JLabel lblProb;
	private JTextField txtHMMProb;
	private JLabel lblRandOdd;
	
	private JScrollPane scrollPane;
	private JPanel panel;
	
	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DnaMLKUserInterface window = new DnaMLKUserInterface(args);
					window.frmDnaMLKControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	protected void ChooseFile(JTextField file) {
		// Construct a new file choose whose default path is the path to this
		// executable, which
		// is returned by System.getProperty("user.dir")

		JFileChooser fileChooser = new JFileChooser(filedir);

		int option = fileChooser.showOpenDialog(frmDnaMLKControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}

	protected void EmpBFToggle(boolean isEmp){
		if (isEmp){
			rdbtnEmpBFYes.setSelected(true);
			rdbtnEmpBFNo.setSelected(false);
			lblBaseFreq.setEnabled(false);
			lblBaseFreqA.setEnabled(false);
			txtBaseFreqA.setEnabled(false);
			lblBaseFreqC.setEnabled(false);
			txtBaseFreqC.setEnabled(false);
			lblBaseFreqG.setEnabled(false);
			txtBaseFreqG.setEnabled(false);
			lblBaseFreqTU.setEnabled(false);
			txtBaseFreqTU.setEnabled(false);
		}
		else{
			rdbtnEmpBFYes.setSelected(false);
			rdbtnEmpBFNo.setSelected(true);
			lblBaseFreq.setEnabled(true);
			lblBaseFreqA.setEnabled(true);
			txtBaseFreqA.setEnabled(true);
			lblBaseFreqC.setEnabled(true);
			txtBaseFreqC.setEnabled(true);
			lblBaseFreqG.setEnabled(true);
			txtBaseFreqG.setEnabled(true);
			lblBaseFreqTU.setEnabled(true);
			txtBaseFreqTU.setEnabled(true);
		}
	}
	
	protected void CatToggle(boolean isOneCat){
		if (isOneCat){
			rdbtnCatYes.setSelected(true);
			rdbtnCatNo.setSelected(false);
			lblNumCat.setEnabled(false);
			cmbxNumCat.setEnabled(false);
			lblcatnum.setEnabled(false);
			ratecatnum.setEnabled(false);
			lblRate.setEnabled(false);
			txtSiteRate.setEnabled(false);
			btnCatFile.setEnabled(false);
			txtCatFile.setEnabled(false);
		}
		else{
			rdbtnCatYes.setSelected(false);
			rdbtnCatNo.setSelected(true);
			lblNumCat.setEnabled(true);
			cmbxNumCat.setEnabled(true);
			lblcatnum.setEnabled(true);
			ratecatnum.setEnabled(true);
			lblRate.setEnabled(true);
			txtSiteRate.setEnabled(true);
			btnCatFile.setEnabled(true);
			txtCatFile.setEnabled(true);
		}
	}
	
	protected void SitesWeightToggle(boolean isSites){
		if (isSites){
			rdbtnSitesWeightYes.setSelected(false);
			rdbtnSitesWeightNo.setSelected(true);
			btnWeightFile.setEnabled(false);
			txtWeightFile.setEnabled(false);
			ExplicitWgts = false;
		}
		else{
			rdbtnSitesWeightYes.setSelected(true);
			rdbtnSitesWeightNo.setSelected(false);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
			ExplicitWgts = true;
		}
	}
	
	protected void GlobalReToggle(boolean isGlobal){
		if(isGlobal){
			rdbtnGlobalReNo.setSelected(true);
			rdbtnGlobalReYes.setSelected(false);
		}
		else{
			rdbtnGlobalReNo.setSelected(false);
			rdbtnGlobalReYes.setSelected(true);
		}
	}
	
	protected void RandOrderToggle(boolean isRand){
		if(isRand){
			rdbtnRandOrderNo.setSelected(false);
			rdbtnRandOrderYes.setSelected(true);
			lblRandSeed.setEnabled(true);
			txtRandSeed.setEnabled(true);
			lblNumberJumble.setEnabled(true);
			txtNumberJumble.setEnabled(true);
			lblRandOdd.setEnabled(true);
		}
		else{
			rdbtnRandOrderNo.setSelected(true);
			rdbtnRandOrderYes.setSelected(false);
			lblRandSeed.setEnabled(false);
			txtRandSeed.setEnabled(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEnabled(false);
			lblRandOdd.setEnabled(false);
		}
	}
	
	protected void MultToggle(boolean isMult){
		if(isMult){
			rdbtnAnalyzeMultDataNo.setSelected(false);
			rdbtnAnalyzeMultDataYes.setSelected(true);
			lblMultData.setEnabled(true);
			rdbtnDataSets.setEnabled(true);
			rdbtnWeights.setEnabled(true);
			lblHowManyData.setEnabled(true);
			txtNumSeqs.setEnabled(true);
			lblInputSeq.setEnabled(true);
			rdbtnInputSeqYes.setEnabled(true);
			rdbtnInputSeqNo.setEnabled(true);
			if (rdbtnWeights.isSelected())
			{
				if(!ExplicitWgts)
				{
					btnWeightFile.setEnabled(true);
					txtWeightFile.setEnabled(true);	
					rdbtnSitesWeightYes.setSelected(true);
					rdbtnSitesWeightNo.setSelected(false);
				}
			}
		}
		else{
			rdbtnAnalyzeMultDataNo.setSelected(true);
			rdbtnAnalyzeMultDataYes.setSelected(false);
			lblMultData.setEnabled(false);
			rdbtnDataSets.setEnabled(false);
			rdbtnWeights.setEnabled(false);
			lblHowManyData.setEnabled(false);
			txtNumSeqs.setEnabled(false);
			lblInputSeq.setEnabled(false);
			rdbtnInputSeqYes.setEnabled(false);
			rdbtnInputSeqNo.setEnabled(false);
			if(ExplicitWgts)
			{
				btnWeightFile.setEnabled(true);
				txtWeightFile.setEnabled(true);	
				rdbtnSitesWeightYes.setSelected(true);
				rdbtnSitesWeightNo.setSelected(false);
			}
			else
			{
				btnWeightFile.setEnabled(false);
				txtWeightFile.setEnabled(false);	
				rdbtnSitesWeightYes.setSelected(false);
				rdbtnSitesWeightNo.setSelected(true);				
			}
	}
	}
	
	protected void InputSeqToggle(boolean isSeq){
		if(isSeq){
			rdbtnInputSeqYes.setSelected(true);
			rdbtnInputSeqNo.setSelected(false);
		}
		else{
			rdbtnInputSeqYes.setSelected(false);
			rdbtnInputSeqNo.setSelected(true);
		}
	}
	
	protected void PrintDataToggle(boolean isData){
		if(isData){
			rdbtnPrintDataYes.setSelected(true);
			rdbtnPrintDataNo.setSelected(false);
		}
		else
		{
			rdbtnPrintDataYes.setSelected(false);
			rdbtnPrintDataNo.setSelected(true);
		}
	}
	
	protected void PrintIndToggle(boolean isInd){
		if(isInd){
			rdbtnPrintIndYes.setSelected(true);
			rdbtnPrintIndNo.setSelected(false);
		}
		else{
			rdbtnPrintIndYes.setSelected(false);
			rdbtnPrintIndNo.setSelected(true);
		}
	}
	
	protected void PrintTreeToggle(boolean isTree){
		if(isTree){
			rdbtnPrintTreeYes.setSelected(true);
			rdbtnPrintTreeNo.setSelected(false);
		}
		else{
			rdbtnPrintTreeYes.setSelected(false);
			rdbtnPrintTreeNo.setSelected(true);
		}
	}
	
	protected void WriteTreeToggle(boolean isWrite){
		if(isWrite){
			rdbtnWriteTreeYes.setSelected(true);
			rdbtnWriteTreeNo.setSelected(false);
			btnOutputTree.setEnabled(true);
			txtOutTree.setEnabled(true);
		}
		else{
			rdbtnWriteTreeYes.setSelected(false);
			rdbtnWriteTreeNo.setSelected(true);
			btnOutputTree.setEnabled(false);
			txtOutTree.setEnabled(false);
		}
	}
	
	protected void RecHypoToggle(boolean isRec){
		if(isRec){
			rdbtnRecHypoNo.setSelected(true);
			rdbtnRecHypoYes.setSelected(false);
		}
		else{
			rdbtnRecHypoNo.setSelected(false);
			rdbtnRecHypoYes.setSelected(true);
		}
	}
	
	protected void DataWeightToggle(boolean isData){
		if(isData){
			rdbtnDataSets.setSelected(true);
			rdbtnWeights.setSelected(false);
			RandOrderToggle(true);
			if(!ExplicitWgts)
			{
				btnWeightFile.setEnabled(false);
				txtWeightFile.setEnabled(false);
				rdbtnSitesWeightYes.setSelected(false);
				rdbtnSitesWeightNo.setSelected(true);
			}
		}
		else{
			rdbtnDataSets.setSelected(false);
			rdbtnWeights.setSelected(true);
			RandOrderToggle(false);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
			rdbtnSitesWeightYes.setSelected(true);
			rdbtnSitesWeightNo.setSelected(false);
		}
	}
	
	protected void IntreeToggle(int selected){
		if(selected == 0){
			txtInputTree.setEnabled(false);
			btnInputTree.setEnabled(false);
			lblUseLengths.setEnabled(false);
			rdbtnUseLengthsYes.setEnabled(false);
			rdbtnUseLengthsNo.setEnabled(false);
		}
		else{
			txtInputTree.setEnabled(true);
			btnInputTree.setEnabled(true);
			lblUseLengths.setEnabled(true);
			rdbtnUseLengthsYes.setEnabled(true);
			rdbtnUseLengthsNo.setEnabled(true);
		}
	}

	protected void DisplaySiteRateValue(int catnum){
		if ((catnum-1) != lastsiteratecat){
			SiteRateValues[lastsiteratecat] = Double.parseDouble(txtSiteRate.getText());
			String siteratestr = Double.toString(SiteRateValues[catnum-1]);
			txtSiteRate.setText(siteratestr);
		}
		lastsiteratecat = catnum-1;
	}
	
	protected void DisplayHMMValues(int catnum){
		if ((catnum-1) != lasthmmcat){
			HMMRateValues[lasthmmcat] = Double.parseDouble(txtHMMRate.getText());
			String hmmratestr = Double.toString(HMMRateValues[catnum-1]);
			txtHMMRate.setText(hmmratestr);
			
			HMMProbValues[lasthmmcat] = Double.parseDouble(txtHMMProb.getText());
			String ratestr = Double.toString(HMMProbValues[catnum-1]);
			txtHMMProb.setText(ratestr);
		}
		lasthmmcat = catnum-1;
	}
	
	protected void RatesiteToggle(int selected){
		if(selected == 0) //Constant rate
		{
			rdbtnAdjCorYes.setEnabled(false);
			rdbtnAdjCorNo.setEnabled(false);
			lblRatesAdjCor.setEnabled(false);
			lblMeanLen.setEnabled(false);
			txtBlockLen.setEnabled(false);
			lblCoeffVarNote.setVisible(false);
			lblInvar.setVisible(false);
			lblHMM.setEnabled(false);
			lblCoeffVar.setEnabled(false);
			txtCoeffVar.setEnabled(false);
			lblNumHMMcats.setEnabled(false);
			cmbxHMMcount.setEnabled(false);
			lblHMMCat.setEnabled(false);
			cmbxHMMCat.setEnabled(false);
			lblHMMRate.setEnabled(false);
			txtHMMRate.setEnabled(false);
			lblFracInvar.setEnabled(false);
			txtFracInvar.setEnabled(false);
			txtHMMProb.setEnabled(false);
			lblProb.setEnabled(false);

		}
		else{
			rdbtnAdjCorYes.setEnabled(true);
			rdbtnAdjCorNo.setEnabled(true);
			lblRatesAdjCor.setEnabled(true);
			lblCoeffVar.setEnabled(true);
			txtCoeffVar.setEnabled(true);
			lblHMM.setEnabled(true);
			lblNumHMMcats.setEnabled(true);
			cmbxHMMcount.setEnabled(true);
			AdjCorToggle(rdbtnAdjCorYes.isSelected());
			if(selected == 1) // Gamma distance rates
			{
				lblCoeffVarNote.setVisible(true);
				lblInvar.setVisible(false);
				lblFracInvar.setEnabled(false);
				txtFracInvar.setEnabled(false);
				txtHMMProb.setEnabled(false);
				lblProb.setEnabled(false);
				lblHMMCat.setEnabled(false);
				cmbxHMMCat.setEnabled(false);
				lblHMMRate.setEnabled(false);
				txtHMMRate.setEnabled(false);
			}
			else if(selected == 2) // Gamma + invariant sites
			{
				lblCoeffVarNote.setVisible(true);
				lblInvar.setText("(including one for invariant sites)");
				lblInvar.setVisible(true);
				lblFracInvar.setEnabled(true);
				txtFracInvar.setEnabled(true);
				txtHMMProb.setEnabled(false);
				lblProb.setEnabled(false);
				lblHMMCat.setEnabled(false);
				cmbxHMMCat.setEnabled(false);
				lblHMMRate.setEnabled(false);
				txtHMMRate.setEnabled(false);
			}
			else  // User-defined HMM of rates
			{
				lblCoeffVarNote.setVisible(false);
				lblInvar.setText("(Probabilities will be scaled to 1.0)");
				lblInvar.setVisible(true);
				lblFracInvar.setEnabled(false);
				txtFracInvar.setEnabled(false);
				txtHMMProb.setEnabled(true);
				lblProb.setEnabled(true);
				lblHMMCat.setEnabled(true);
				cmbxHMMCat.setEnabled(true);
				lblHMMRate.setEnabled(true);
				txtHMMRate.setEnabled(true);
			}
		}
	}
	
	protected void AdjCorToggle(boolean isYes)
	{
		if(isYes)
		{
			rdbtnAdjCorYes.setSelected(true);
			rdbtnAdjCorNo.setSelected(false);
			lblMeanLen.setEnabled(true);
			txtBlockLen.setEnabled(true);			
		}
		else
		{
			rdbtnAdjCorYes.setSelected(false);
			rdbtnAdjCorNo.setSelected(true);
			lblMeanLen.setEnabled(false);
			txtBlockLen.setEnabled(false);			
		}	
	}

	protected void DotDiffToggle(boolean isUseDot) {
		if (isUseDot) {
			rdbtnDotDiffYes.setSelected(true);
			rdbtnDotDiffNo.setSelected(false);
		} else {
			rdbtnDotDiffYes.setSelected(false);
			rdbtnDotDiffNo.setSelected(true);
		}
	}

	protected void UseLengthsToggle(boolean isUseLength) {
		if (isUseLength) {
			rdbtnUseLengthsYes.setSelected(true);
			rdbtnUseLengthsNo.setSelected(false);
		} else {
			rdbtnUseLengthsYes.setSelected(false);
			rdbtnUseLengthsNo.setSelected(true);
		}
	}

	/**
	 * Create the application.
	 */
	public DnaMLKUserInterface(String[] args) {
		phylipCall = false;
		if (args.length > 0)
		{
			if (args[0].contains("Phylip"))
			{
				phylipCall = true;
			}
		}
		initialize();
	}
	
	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		
		SiteRateValues = new double[9];
		SiteRateValues[0] = 1.0;
		SiteRateValues[1] = 1.0;
		SiteRateValues[2] = 1.0;
		SiteRateValues[3] = 1.0;
		SiteRateValues[4] = 1.0;
		SiteRateValues[5] = 1.0;
		SiteRateValues[6] = 1.0;
		SiteRateValues[7] = 1.0;
		SiteRateValues[8] = 1.0;
		
		lastsiteratecat = 0;
		
		HMMRateValues = new double[9];
		HMMRateValues[0] = 1.0;
		HMMRateValues[1] = 1.0;
		HMMRateValues[2] = 1.0;
		HMMRateValues[3] = 1.0;
		HMMRateValues[4] = 1.0;
		HMMRateValues[5] = 1.0;
		HMMRateValues[6] = 1.0;
		HMMRateValues[7] = 1.0;
		HMMRateValues[8] = 1.0;
		
		lasthmmcat = 0;
		
		HMMProbValues = new double[9];
		HMMProbValues[0] = 1.0;
		HMMProbValues[1] = 1.0;
		HMMProbValues[2] = 1.0;
		HMMProbValues[3] = 1.0;
		HMMProbValues[4] = 1.0;
		HMMProbValues[5] = 1.0;
		HMMProbValues[6] = 1.0;
		HMMProbValues[7] = 1.0;
		HMMProbValues[8] = 1.0;
		
		filedir = System.getProperty("user.dir");
		ExplicitWgts = false;

		frmDnaMLKControls = new JFrame();
		frmDnaMLKControls.setBackground(new Color(204, 255, 255));
		frmDnaMLKControls.setTitle("Dnamlk");
		frmDnaMLKControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDnaMLKControls.setBounds(100, 100, 1200, 750);
		frmDnaMLKControls.setPreferredSize(new Dimension(frmDnaMLKControls.getBounds().width, frmDnaMLKControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmDnaMLKControls.getPreferredSize());
		frmDnaMLKControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmDnaMLKControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow][pref!,grow][pref!,grow]", "[][][][]"));
		
		btnInputFile = new JButton("Input File");
		btnInputFile.setBounds(21, 10, 141, 23);
		btnInputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputFile);
			}
		});
		panel.add(btnInputFile, "cell 0 0,growx");
		
		txtInputFile = new JTextField();
		txtInputFile.setBounds(168, 11, 852, 20);
		txtInputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInputFile.setText("infile");
		panel.add(txtInputFile, "cell 1 0 5 1,growx");
		
		btnInputTree = new JButton("Input Tree");
		btnInputTree.setBounds(21, 35, 141, 23);
		btnInputTree.setEnabled(false);
		btnInputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputTree);
			}
		});
		btnInputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnInputTree, "cell 0 1,growx");
		
		txtInputTree = new JTextField();
		txtInputTree.setBounds(168, 36, 852, 20);
		txtInputTree.setEnabled(false);
		txtInputTree.setText("intree");
		txtInputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInputTree, "cell 1 1 5 1,growx");
		
		btnWeightFile = new JButton("Weights File");
		btnWeightFile.setBounds(21, 60, 141, 23);
		btnWeightFile.setEnabled(false);
		btnWeightFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightFile);
			}
		});
		btnWeightFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnWeightFile, "cell 0 2,growx");
	
		txtWeightFile = new JTextField();
		txtWeightFile.setBounds(168, 61, 852, 20);
		txtWeightFile.setEnabled(false);
		txtWeightFile.setText("weightfile");
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtWeightFile, "cell 1 2 5 1,growx");
		
		btnCatFile = new JButton("Categories File");
		btnCatFile.setBounds(21, 85, 141, 23);
		btnCatFile.setEnabled(false);
		btnCatFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtCatFile);
			}
		});
		btnCatFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnCatFile, "cell 0 3,growx");
		
		txtCatFile = new JTextField();
		txtCatFile.setBounds(168, 86, 852, 20);
		txtCatFile.setEnabled(false);
		txtCatFile.setText("catfile");
		txtCatFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtCatFile, "cell 1 3 5 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.setBounds(21, 110, 141, 23);
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		panel.add(btnOutputFile, "cell 0 4,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setBounds(168, 111, 852, 20);
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 4 5 1,growx");

		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.setBounds(21, 135, 141, 23);
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutTree);
			}
		});
		panel.add(btnOutputTree, "cell 0 5,growx");
		
		txtOutTree = new JTextField();
		txtOutTree.setBounds(168, 136, 852, 20);
		txtOutTree.setText("outtree");
		txtOutTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutTree, "cell 1 5 5 1,growx");
	
		lblSearchBest = new JLabel("Search for best tree:");
		lblSearchBest.setBounds(114, 175, 141, 14);
		lblSearchBest.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSearchBest.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSearchBest, "flowx,cell 0 6 2 1,alignx right");
		
		TreeSearchMethod = new JComboBox();
		TreeSearchMethod.setBounds(267, 171, 264, 23);
		TreeSearchMethod.setModel(new DefaultComboBoxModel(new String[] {"Yes", "No, use user trees in input file"}));
		TreeSearchMethod.setSelectedIndex(0);
		TreeSearchMethod.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(TreeSearchMethod.getSelectedIndex());
			}
		});
		panel.add(TreeSearchMethod, "cell 2 6 2 1,growx");
		
		lblUseLengths = new JLabel("Use lengths from user trees:");
		lblUseLengths.setBounds(36, 200, 219, 14);
		lblUseLengths.setEnabled(false);
		lblUseLengths.setHorizontalAlignment(SwingConstants.RIGHT);
		lblUseLengths.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseLengths, "flowx,cell 0 7 2 1,alignx right");
		
		rdbtnUseLengthsYes = new JRadioButton("Yes");
		rdbtnUseLengthsYes.setBounds(267, 196, 56, 23);
		rdbtnUseLengthsYes.setEnabled(false);
		rdbtnUseLengthsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseLengthsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UseLengthsToggle(true);
			}
		});
		rdbtnUseLengthsYes.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnUseLengthsYes, "cell 2 7");
		
		rdbtnUseLengthsNo = new JRadioButton("No");
		rdbtnUseLengthsNo.setBounds(328, 196, 55, 23);
		rdbtnUseLengthsNo.setSelected(true);
		rdbtnUseLengthsNo.setEnabled(false);
		rdbtnUseLengthsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseLengthsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UseLengthsToggle(false);
			}
		});
		rdbtnUseLengthsNo.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnUseLengthsNo, "cell 2 7");
		
		lblTTratio = new JLabel("Transition/transversion ratio:");
		lblTTratio.setBounds(66, 225, 190, 14);
		lblTTratio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblTTratio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTTratio, "flowx,cell 0 8 2 1,alignx right");
		
		txtTTratio = new JTextField();
		txtTTratio.setBounds(268, 222, 59, 20);
		txtTTratio.setFont(new Font("Arial", Font.PLAIN, 13));
		txtTTratio.setText("2.00");
		txtTTratio.setColumns(5);
		panel.add(txtTTratio, "cell 2 8");
		
		lblEmpBF = new JLabel("Use empirical base frequencies:");
		lblEmpBF.setBounds(37, 250, 219, 14);
		lblEmpBF.setHorizontalAlignment(SwingConstants.RIGHT);
		lblEmpBF.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblEmpBF, "flowx,cell 0 9 2 1,alignx right");
		
		rdbtnEmpBFYes = new JRadioButton("Yes");
		rdbtnEmpBFYes.setBounds(268, 246, 56, 23);
		rdbtnEmpBFYes.setSelected(true);
		rdbtnEmpBFYes.setBackground(new Color(204, 255, 255));
		rdbtnEmpBFYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnEmpBFYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EmpBFToggle(true);
			}
		});
		panel.add(rdbtnEmpBFYes, "cell 2 9");
		
		rdbtnEmpBFNo = new JRadioButton("No");
		rdbtnEmpBFNo.setBounds(326, 246, 55, 23);
		rdbtnEmpBFNo.setBackground(new Color(204, 255, 255));
		rdbtnEmpBFNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnEmpBFNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EmpBFToggle(false);
			}
		});
		panel.add(rdbtnEmpBFNo, "cell 2 9");
		
		lblBaseFreq = new JLabel("Base frequencies:");
		lblBaseFreq.setBounds(105, 275, 151, 14);
		lblBaseFreq.setEnabled(false);
		lblBaseFreq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblBaseFreq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreq, "flowx,cell 0 10 2 1,alignx right");
		
		lblBaseFreqA = new JLabel("A:");
		lblBaseFreqA.setBounds(268, 275, 23, 14);
		lblBaseFreqA.setEnabled(false);
		lblBaseFreqA.setBackground(new Color(153, 255, 255));
		lblBaseFreqA.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqA, "cell 2 10");
		
		txtBaseFreqA = new JTextField();
		txtBaseFreqA.setBounds(281, 272, 41, 20);
		txtBaseFreqA.setEnabled(false);
		txtBaseFreqA.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqA.setText("0.25");
		txtBaseFreqA.setColumns(6);
		panel.add(txtBaseFreqA, "cell 2 10");
		
		lblBaseFreqC = new JLabel("   C:");
		lblBaseFreqC.setBounds(343, 275, 23, 14);
		lblBaseFreqC.setEnabled(false);
		lblBaseFreqC.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqC, "cell 3 10");
		
		txtBaseFreqC = new JTextField();
		txtBaseFreqC.setBounds(358, 272, 41, 20);
		txtBaseFreqC.setEnabled(false);
		txtBaseFreqC.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqC.setText("0.25");
		txtBaseFreqC.setColumns(6);
		panel.add(txtBaseFreqC, "cell 3 10");
		
		lblBaseFreqG = new JLabel("G:");
		lblBaseFreqG.setBounds(268, 304, 23, 14);
		lblBaseFreqG.setEnabled(false);
		lblBaseFreqG.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqG, "cell 2 11");
		
		txtBaseFreqG = new JTextField();
		txtBaseFreqG.setBounds(282, 301, 41, 20);
		txtBaseFreqG.setEnabled(false);
		txtBaseFreqG.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqG.setText("0.25");
		txtBaseFreqG.setColumns(6);
		panel.add(txtBaseFreqG, "cell 2 11");
		
		lblBaseFreqTU = new JLabel("T/U:");
		lblBaseFreqTU.setBounds(338, 304, 23, 14);
		lblBaseFreqTU.setEnabled(false);
		lblBaseFreqTU.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqTU, "cell 3 11");
		
		txtBaseFreqTU = new JTextField();
		txtBaseFreqTU.setBounds(362, 301, 41, 20);
		txtBaseFreqTU.setEnabled(false);
		txtBaseFreqTU.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqTU.setText("0.25");
		txtBaseFreqTU.setColumns(6);
		panel.add(txtBaseFreqTU, "cell 3 11");
		
		lblCatSites = new JLabel("One category of substitution rates:");
		lblCatSites.setBounds(37, 328, 220, 19);
		lblCatSites.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCatSites.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCatSites, "flowx,cell 0 12 2 1,alignx right");
		
		rdbtnCatYes = new JRadioButton("Yes");
		rdbtnCatYes.setBounds(266, 325, 56, 23);
		rdbtnCatYes.setSelected(true);
		rdbtnCatYes.setBackground(new Color(204, 255, 255));
		rdbtnCatYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCatYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CatToggle(true);
			}
		});
		panel.add(rdbtnCatYes, "cell 2 12");
		
		rdbtnCatNo = new JRadioButton("No");
		rdbtnCatNo.setBounds(324, 325, 55, 23);
		rdbtnCatNo.setBackground(new Color(204, 255, 255));
		rdbtnCatNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCatNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CatToggle(false);
			}
		});
		panel.add(rdbtnCatNo, "cell 2 12");
		
		lblNumCat = new JLabel("Number of substitution categories:");
		lblNumCat.setBounds(21, 351, 235, 17);
		lblNumCat.setEnabled(false);
		lblNumCat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumCat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNumCat, "flowx,cell 0 13 2 1,alignx right");
		
		cmbxNumCat = new JComboBox();
		cmbxNumCat.setBounds(266, 349, 59, 23);
		cmbxNumCat.setEnabled(false);
		cmbxNumCat.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxNumCat.setSelectedIndex(1);
		cmbxNumCat.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxNumCat, "cell 2 13");
		
		lblcatnum = new JLabel("Category:");
		lblcatnum.setBounds(332, 350, 61, 20);
		lblcatnum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblcatnum.setFont(new Font("Arial", Font.BOLD, 13));
		lblcatnum.setEnabled(false);
		panel.add(lblcatnum, "cell 2 13");
		
		ratecatnum = new JComboBox();
		ratecatnum.setBounds(396, 349, 61, 23);
		ratecatnum.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		ratecatnum.setFont(new Font("Arial", Font.PLAIN, 13));
		ratecatnum.setEnabled(false);
		ratecatnum.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplaySiteRateValue(Integer.parseInt(ratecatnum.getSelectedItem().toString()));
			}
		});
		panel.add(ratecatnum, "cell 3 13");
		
		lblRate = new JLabel("Rate:");
		lblRate.setBounds(456, 350, 35, 20);
		lblRate.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRate.setFont(new Font("Arial", Font.BOLD, 13));
		lblRate.setEnabled(false);
		panel.add(lblRate, "cell 3 13");
		
		txtSiteRate = new JTextField();
		txtSiteRate.setBounds(495, 350, 51, 20);
		txtSiteRate.setText("1.0");
		txtSiteRate.setHorizontalAlignment(SwingConstants.CENTER);
		txtSiteRate.setFont(new Font("Arial", Font.PLAIN, 13));
		txtSiteRate.setEnabled(false);
		txtSiteRate.setColumns(5);
		panel.add(txtSiteRate, "cell 3 13");
	
		lblRateSite = new JLabel("Rate variation among sites:");
		lblRateSite.setBounds(64, 377, 190, 17);
		lblRateSite.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRateSite.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRateSite, "flowx,cell 0 14 2 1,alignx right");
		
		cmbxRateSite = new JComboBox();
		cmbxRateSite.setBounds(266, 374, 282, 23);
		cmbxRateSite.setModel(new DefaultComboBoxModel(new String[] {"Constant rate", "Gamma distributed rates", "Gamma + invariant sites", "User-defined HMM of rates"}));
		cmbxRateSite.setSelectedIndex(0);
		cmbxRateSite.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxRateSite.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				RatesiteToggle(cmbxRateSite.getSelectedIndex());
			}
		});
		panel.add(cmbxRateSite, "cell 2 14 2 1,growx");
		
		lblRatesAdjCor = new JLabel("Rates on adjacent sites correlated:");
		lblRatesAdjCor.setBounds(28, 402, 226, 16);
		lblRatesAdjCor.setEnabled(false);
		lblRatesAdjCor.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRatesAdjCor.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRatesAdjCor, "flowx,cell 0 15 2 1,alignx right");
		
		rdbtnAdjCorYes = new JRadioButton("Yes");
		rdbtnAdjCorYes.setBounds(266, 399, 56, 23);
		rdbtnAdjCorYes.setEnabled(false);
		rdbtnAdjCorYes.setSelected(false);
		rdbtnAdjCorYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAdjCorYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AdjCorToggle(true);
			}
		});
		rdbtnAdjCorYes.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnAdjCorYes, "cell 2 15");
		
		rdbtnAdjCorNo = new JRadioButton("No");
		rdbtnAdjCorNo.setBounds(324, 399, 64, 23);
		rdbtnAdjCorNo.setSelected(true);
		rdbtnAdjCorNo.setEnabled(false);
		rdbtnAdjCorNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAdjCorNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AdjCorToggle(false);
			}
		});
		rdbtnAdjCorNo.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnAdjCorNo, "cell 2 15");
		
		lblMeanLen = new JLabel("Mean block length:");
		lblMeanLen.setBounds(372, 400, 128, 20);
		lblMeanLen.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMeanLen.setFont(new Font("Arial", Font.BOLD, 13));
		lblMeanLen.setEnabled(false);
		panel.add(lblMeanLen, "cell 3 15");
		
		txtBlockLen = new JTextField();
		txtBlockLen.setBounds(512, 400, 43, 20);
		txtBlockLen.setText("1");
		txtBlockLen.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBlockLen.setEnabled(false);
		txtBlockLen.setColumns(5);
		panel.add(txtBlockLen, "cell 3 15");
		
		lblCoeffVar = new JLabel("Coefficient of variation:");
		lblCoeffVar.setBounds(94, 427, 159, 14);
		lblCoeffVar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCoeffVar.setFont(new Font("Arial", Font.BOLD, 13));
		lblCoeffVar.setEnabled(false);
		panel.add(lblCoeffVar, "flowx,cell 0 16 2 1,alignx right");

		txtCoeffVar = new JTextField();
		txtCoeffVar.setBounds(266, 424, 75, 20);
		txtCoeffVar.setEnabled(false);
		txtCoeffVar.setText("1.0");
		txtCoeffVar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCoeffVar.setColumns(6);
		panel.add(txtCoeffVar, "cell 2 16");
		
		lblCoeffVarNote = new JLabel("(for gamma dist = 1/\u221Aalpha)");
		lblCoeffVarNote.setBounds(358, 424, 262, 14);
		lblCoeffVarNote.setHorizontalAlignment(SwingConstants.LEFT);
		lblCoeffVarNote.setFont(new Font("Arial", Font.BOLD, 13));
		lblCoeffVarNote.setVisible(false);
		panel.add(lblCoeffVarNote, "cell 3 16");
		
		lblHMM = new JLabel("Rates for HMM:");
		lblHMM.setBounds(95, 452, 159, 14);
		lblHMM.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMM.setFont(new Font("Arial", Font.BOLD, 13));
		lblHMM.setEnabled(false);
		panel.add(lblHMM, "flowx,cell 0 17 2 1,alignx right");
		
		lblInvar = new JLabel("");
		lblInvar.setBounds(268, 422, 328, 14);
		lblInvar.setHorizontalAlignment(SwingConstants.LEFT);
		lblInvar.setFont(new Font("Arial", Font.BOLD, 13));
		lblInvar.setVisible(false);
		panel.add(lblInvar, "cell 2 17");
		
		lblNumHMMcats = new JLabel("Number of rate categories:");
		lblNumHMMcats.setBounds(78, 477, 176, 17);
		lblNumHMMcats.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumHMMcats.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumHMMcats.setEnabled(false);
		panel.add(lblNumHMMcats, "flowx,cell 0 18 2 1,alignx right");

		cmbxHMMcount = new JComboBox();
		cmbxHMMcount.setBounds(266, 474, 59, 23);
		cmbxHMMcount.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxHMMcount.setSelectedIndex(0);
		cmbxHMMcount.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxHMMcount.setEnabled(false);
		panel.add(cmbxHMMcount, "cell 2 18,growx");
		
		lblHMMCat = new JLabel("Category:");
		lblHMMCat.setBounds(332, 475, 61, 20);
		lblHMMCat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMMCat.setFont(new Font("Arial", Font.BOLD, 13));
		lblHMMCat.setEnabled(false);
		panel.add(lblHMMCat, "cell 3 18");
		
		cmbxHMMCat = new JComboBox();
		cmbxHMMCat.setBounds(396, 474, 61, 23);
		cmbxHMMCat.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxHMMCat.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxHMMCat.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplayHMMValues(Integer.parseInt(cmbxHMMCat.getSelectedItem().toString()));
			}
		});
		cmbxHMMCat.setEnabled(false);
		panel.add(cmbxHMMCat, "cell 3 18,growx");
		
		lblHMMRate = new JLabel("Rate:");
		lblHMMRate.setBounds(279, 499, 35, 20);
		lblHMMRate.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMMRate.setFont(new Font("Arial", Font.BOLD, 13));
		lblHMMRate.setEnabled(false);
		panel.add(lblHMMRate, "cell 2 19");
		
		txtHMMRate = new JTextField();
		txtHMMRate.setBounds(318, 499, 51, 20);
		txtHMMRate.setText("1.0");
		txtHMMRate.setHorizontalAlignment(SwingConstants.CENTER);
		txtHMMRate.setFont(new Font("Arial", Font.PLAIN, 13));
		txtHMMRate.setEnabled(false);
		txtHMMRate.setColumns(6);
		panel.add(txtHMMRate, "cell 2 19");
		
		lblProb = new JLabel("Probablilty:");
		lblProb.setBounds(371, 499, 75, 20);
		lblProb.setHorizontalAlignment(SwingConstants.RIGHT);
		lblProb.setFont(new Font("Arial", Font.BOLD, 13));
		lblProb.setEnabled(false);
		panel.add(lblProb, "cell 3 19");
		
		txtHMMProb = new JTextField();
		txtHMMProb.setBounds(450, 499, 51, 20);
		txtHMMProb.setText("1.0");
		txtHMMProb.setHorizontalAlignment(SwingConstants.CENTER);
		txtHMMProb.setFont(new Font("Arial", Font.PLAIN, 13));
		txtHMMProb.setEnabled(false);
		txtHMMProb.setColumns(6);
		panel.add(txtHMMProb, "cell 3 19");
	
		lblFracInvar = new JLabel("Fraction of sites invariant:");
		lblFracInvar.setBounds(64, 525, 189, 14);
		lblFracInvar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblFracInvar.setFont(new Font("Arial", Font.BOLD, 13));
		lblFracInvar.setEnabled(false);
		panel.add(lblFracInvar, "flowx,cell 0 20 2 1,alignx right");
	
		txtFracInvar = new JTextField();
		txtFracInvar.setBounds(266, 522, 75, 20);
		txtFracInvar.setText("1.0");
		txtFracInvar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtFracInvar.setEnabled(false);
		txtFracInvar.setColumns(6);
		panel.add(txtFracInvar, "cell 2 20");

		// **** column 2 ****//	
		lblGlobalRe = new JLabel("Global rearrangements:");
		lblGlobalRe.setBounds(647, 169, 159, 20);
		lblGlobalRe.setHorizontalAlignment(SwingConstants.RIGHT);
		lblGlobalRe.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGlobalRe, "cell 4 6,alignx right");
		
		rdbtnGlobalReYes = new JRadioButton("Yes");
		rdbtnGlobalReYes.setBounds(819, 168, 56, 23);
		rdbtnGlobalReYes.setBackground(new Color(204, 255, 255));
		rdbtnGlobalReYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnGlobalReYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalReToggle(false);
			}
		});
		rdbtnGlobalReYes.setSelected(false);
		panel.add(rdbtnGlobalReYes, "cell 5 6");
		
		rdbtnGlobalReNo = new JRadioButton("No");
		rdbtnGlobalReNo.setBounds(877, 168, 55, 23);
		rdbtnGlobalReNo.setBackground(new Color(204, 255, 255));
		rdbtnGlobalReNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnGlobalReNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalReToggle(true);
			}
		});
		rdbtnGlobalReNo.setSelected(true);
		panel.add(rdbtnGlobalReNo, "cell 5 6");

		lblSitesWeight = new JLabel("Sites weighted:");
		lblSitesWeight.setBounds(685, 194, 121, 16);
		lblSitesWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSitesWeight.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSitesWeight, "cell 4 7,alignx right");
		
		rdbtnSitesWeightYes = new JRadioButton("Yes");
		rdbtnSitesWeightYes.setBounds(819, 191, 56, 23);
		rdbtnSitesWeightYes.setBackground(new Color(204, 255, 255));
		rdbtnSitesWeightYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesWeightYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SitesWeightToggle(false);
			}
		});
		rdbtnSitesWeightYes.setSelected(false);
		panel.add(rdbtnSitesWeightYes, "cell 5 7");
		
		rdbtnSitesWeightNo = new JRadioButton("No");
		rdbtnSitesWeightNo.setBounds(877, 191, 64, 23);
		rdbtnSitesWeightNo.setBackground(new Color(204, 255, 255));
		rdbtnSitesWeightNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesWeightNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SitesWeightToggle(true);
			}
		});
		rdbtnSitesWeightNo.setSelected(true);
		panel.add(rdbtnSitesWeightNo, "cell 5 7");
		
		lblRandOrder = new JLabel("Randomize input order of sequences:");
		lblRandOrder.setBounds(560, 219, 246, 17);
		lblRandOrder.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandOrder.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOrder, "cell 4 8,alignx right");
		
		rdbtnRandOrderYes = new JRadioButton("Yes");
		rdbtnRandOrderYes.setBounds(819, 216, 56, 23);
		rdbtnRandOrderYes.setBackground(new Color(204, 255, 255));
		rdbtnRandOrderYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRandOrderYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(true);
			}
		});
		rdbtnRandOrderYes.setSelected(false);
		panel.add(rdbtnRandOrderYes, "cell 5 8");
		
		rdbtnRandOrderNo = new JRadioButton("No, use input order");
		rdbtnRandOrderNo.setBounds(877, 216, 159, 23);
		rdbtnRandOrderNo.setBackground(new Color(204, 255, 255));
		rdbtnRandOrderNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRandOrderNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(false);
			}
		});
		rdbtnRandOrderNo.setSelected(true);
		panel.add(rdbtnRandOrderNo, "cell 5 8");
		
		lblRandSeed = new JLabel("Random number seed:");
		lblRandSeed.setBounds(652, 244, 154, 17);
		lblRandSeed.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandSeed.setFont(new Font("Arial", Font.BOLD, 13));
		lblRandSeed.setEnabled(false);
		panel.add(lblRandSeed, "cell 4 9,alignx right");
		
		txtRandSeed = new JTextField();
		txtRandSeed.setBounds(819, 242, 86, 20);
		txtRandSeed.setText("1");
		txtRandSeed.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRandSeed.setEnabled(false);
		txtRandSeed.setColumns(6);	
		panel.add(txtRandSeed, "cell 5 9");
		
		lblRandOdd = new JLabel("(must be odd)");
		lblRandOdd.setBounds(911, 245, 109, 14);
		lblRandOdd.setEnabled(false);
		lblRandOdd.setHorizontalAlignment(SwingConstants.LEFT);
		lblRandOdd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOdd, "cell 5 9");

		lblNumberJumble = new JLabel("Number of times to jumble:");
		lblNumberJumble.setBounds(586, 269, 220, 23);
		lblNumberJumble.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberJumble.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumberJumble.setEnabled(false);
		panel.add(lblNumberJumble, "cell 4 10,alignx right");

		txtNumberJumble = new JTextField();
		txtNumberJumble.setBounds(819, 270, 86, 20);
		txtNumberJumble.setText("1");
		txtNumberJumble.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberJumble.setEnabled(false);
		txtNumberJumble.setColumns(6);
		panel.add(txtNumberJumble, "cell 5 10");
	
		lblAnalyzeMultData = new JLabel("Analyze multiple data sets:");
		lblAnalyzeMultData.setBounds(619, 294, 187, 14);
		lblAnalyzeMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAnalyzeMultData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAnalyzeMultData, "cell 4 11,alignx right");
		
		rdbtnAnalyzeMultDataYes = new JRadioButton("Yes");
		rdbtnAnalyzeMultDataYes.setBounds(819, 290, 56, 23);
		rdbtnAnalyzeMultDataYes.setBackground(new Color(204, 255, 255));
		rdbtnAnalyzeMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAnalyzeMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		rdbtnAnalyzeMultDataYes.setSelected(false);
		panel.add(rdbtnAnalyzeMultDataYes, "cell 5 11");
		
		rdbtnAnalyzeMultDataNo = new JRadioButton("No");
		rdbtnAnalyzeMultDataNo.setBounds(877, 290, 50, 23);
		rdbtnAnalyzeMultDataNo.setBackground(new Color(204, 255, 255));
		rdbtnAnalyzeMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAnalyzeMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		rdbtnAnalyzeMultDataNo.setSelected(true);
		panel.add(rdbtnAnalyzeMultDataNo, "cell 5 11");
		
		lblMultData = new JLabel("Multiple data sets or multiple weights:");
		lblMultData.setBounds(565, 319, 241, 19);
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultData.setEnabled(false);
		panel.add(lblMultData, "cell 4 12,alignx right");
		
		rdbtnDataSets = new JRadioButton("Data sets");
		rdbtnDataSets.setBounds(819, 317, 92, 23);
		rdbtnDataSets.setBackground(new Color(204, 255, 255));
		rdbtnDataSets.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(true);
			}
		});
		rdbtnDataSets.setSelected(true);
		rdbtnDataSets.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDataSets.setEnabled(false);
		panel.add(rdbtnDataSets, "cell 5 12");
		
		rdbtnWeights = new JRadioButton("Weights");
		rdbtnWeights.setBounds(926, 317, 92, 23);
		rdbtnWeights.setBackground(new Color(204, 255, 255));
		rdbtnWeights.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(false);
			}
		});
		
		rdbtnWeights.setSelected(false);
		rdbtnWeights.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWeights.setEnabled(false);
		panel.add(rdbtnWeights, "cell 5 12");

		lblHowManyData = new JLabel("Number:");
		lblHowManyData.setBounds(819, 349, 61, 20);
		lblHowManyData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHowManyData.setFont(new Font("Arial", Font.BOLD, 13));
		lblHowManyData.setEnabled(false);
		panel.add(lblHowManyData, "cell 5 13");
		
		txtNumSeqs = new JTextField();
		txtNumSeqs.setBounds(882, 349, 43, 20);
		txtNumSeqs.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumSeqs.setText("1");
		txtNumSeqs.setEnabled(false);
		txtNumSeqs.setColumns(6);
		panel.add(txtNumSeqs, "cell 5 13");
		
		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setBounds(626, 385, 187, 19);
		lblInputSeq.setEnabled(false);
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblInputSeq, "cell 4 14,alignx right");
		
		rdbtnInputSeqYes = new JRadioButton("Interleaved");
		rdbtnInputSeqYes.setBounds(826, 383, 113, 23);
		rdbtnInputSeqYes.setEnabled(false);
		rdbtnInputSeqYes.setBackground(new Color(204, 255, 255));
		rdbtnInputSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		rdbtnInputSeqYes.setSelected(true);
		panel.add(rdbtnInputSeqYes, "cell 5 14");
		
		rdbtnInputSeqNo = new JRadioButton("Sequential");
		rdbtnInputSeqNo.setBounds(943, 383, 121, 23);
		rdbtnInputSeqNo.setEnabled(false);
		rdbtnInputSeqNo.setBackground(new Color(204, 255, 255));
		rdbtnInputSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		rdbtnInputSeqNo.setSelected(false);
		panel.add(rdbtnInputSeqNo, "cell 5 14");
		
		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setBounds(596, 418, 219, 17);
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintData, "cell 4 15,alignx right");
		
		rdbtnPrintDataYes = new JRadioButton("Yes");
		rdbtnPrintDataYes.setBounds(828, 415, 56, 23);
		rdbtnPrintDataYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		rdbtnPrintDataYes.setSelected(false);
		panel.add(rdbtnPrintDataYes, "cell 5 15");
		
		rdbtnPrintDataNo = new JRadioButton("No");
		rdbtnPrintDataNo.setBounds(886, 415, 50, 23);
		rdbtnPrintDataNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rdbtnPrintDataNo.setSelected(true);
		panel.add(rdbtnPrintDataNo, "cell 5 15");
		
		lblPrintTree = new JLabel("Print out tree:");
		lblPrintTree.setBounds(693, 443, 122, 14);
		lblPrintTree.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintTree, "cell 4 16,alignx right");
		
		rdbtnPrintTreeYes = new JRadioButton("Yes");
		rdbtnPrintTreeYes.setBounds(828, 439, 56, 23);
		rdbtnPrintTreeYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(true);
			}
		});
		rdbtnPrintTreeYes.setSelected(true);
		panel.add(rdbtnPrintTreeYes, "cell 5 16");
		
		rdbtnPrintTreeNo = new JRadioButton("No");
		rdbtnPrintTreeNo.setBounds(886, 439, 50, 23);
		rdbtnPrintTreeNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(false);
			}
		});
		rdbtnPrintTreeNo.setSelected(false);
		panel.add(rdbtnPrintTreeNo, "cell 5 16");
		
		lblWriteTree = new JLabel("Write out trees onto tree file:");
		lblWriteTree.setBounds(625, 468, 190, 14);
		lblWriteTree.setHorizontalAlignment(SwingConstants.RIGHT);
		lblWriteTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblWriteTree, "cell 4 17,alignx right");
		
		rdbtnWriteTreeYes = new JRadioButton("Yes");
		rdbtnWriteTreeYes.setBounds(828, 464, 56, 23);
		rdbtnWriteTreeYes.setBackground(new Color(204, 255, 255));
		rdbtnWriteTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(true);
			}
		});
		rdbtnWriteTreeYes.setSelected(true);
		panel.add(rdbtnWriteTreeYes, "cell 5 17");
		
		rdbtnWriteTreeNo = new JRadioButton("No");
		rdbtnWriteTreeNo.setBounds(886, 464, 55, 23);
		rdbtnWriteTreeNo.setBackground(new Color(204, 255, 255));
		rdbtnWriteTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(false);
			}
		});
		rdbtnWriteTreeNo.setSelected(false);
		panel.add(rdbtnWriteTreeNo, "cell 5 17");

		lblDotDiff = new JLabel("Use dot-differencing to display them:");
		lblDotDiff.setBounds(572, 493, 246, 19);
		lblDotDiff.setFont(new Font("Arial", Font.BOLD, 13));
		lblDotDiff.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblDotDiff, "cell 4 18,alignx right");

		rdbtnDotDiffYes = new JRadioButton("Yes");
		rdbtnDotDiffYes.setBounds(828, 491, 56, 23);
		rdbtnDotDiffYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(true);
			}
		});
		rdbtnDotDiffYes.setSelected(true);
		panel.add(rdbtnDotDiffYes, "cell 5 18");

		rdbtnDotDiffNo = new JRadioButton("No");
		rdbtnDotDiffNo.setBounds(886, 491, 55, 23);
		rdbtnDotDiffNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(false);
			}
		});
		panel.add(rdbtnDotDiffNo, "cell 5 18");
		
		lblRecHypo = new JLabel("Reconstruct hypothetical sequences:");
		lblRecHypo.setBounds(581, 518, 234, 14);
		lblRecHypo.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRecHypo.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRecHypo, "cell 4 19,alignx right");
		
		rdbtnRecHypoYes = new JRadioButton("Yes");
		rdbtnRecHypoYes.setBounds(828, 514, 56, 23);
		rdbtnRecHypoYes.setBackground(new Color(204, 255, 255));
		rdbtnRecHypoYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRecHypoYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RecHypoToggle(false);
			}
		});
		rdbtnRecHypoYes.setSelected(false);
		panel.add(rdbtnRecHypoYes, "cell 5 19");
		
		rdbtnRecHypoNo = new JRadioButton("No");
		rdbtnRecHypoNo.setBounds(886, 514, 55, 23);
		rdbtnRecHypoNo.setBackground(new Color(204, 255, 255));
		rdbtnRecHypoNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRecHypoNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RecHypoToggle(true);
			}
		});
		rdbtnRecHypoNo.setSelected(true);
		panel.add(rdbtnRecHypoNo, "cell 5 19");
		
		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setBounds(581, 543, 234, 18);
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInd, "cell 4 20,alignx right");
		
		rdbtnPrintIndYes = new JRadioButton("Yes");
		rdbtnPrintIndYes.setBounds(828, 541, 56, 23);
		rdbtnPrintIndYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rdbtnPrintIndYes.setSelected(true);
		panel.add(rdbtnPrintIndYes, "cell 5 20");
		
		rdbtnPrintIndNo = new JRadioButton("No");
		rdbtnPrintIndNo.setBounds(886, 541, 55, 23);
		rdbtnPrintIndNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rdbtnPrintIndNo.setSelected(false);
		panel.add(rdbtnPrintIndNo, "cell 5 20");
		
		btnExecute = new JButton("Execute");
		btnExecute.setBounds(860, 576, 89, 23);
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				// catch last array data entries
				SiteRateValues[Integer.parseInt(ratecatnum.getSelectedItem().toString())-1] = Double.parseDouble(txtSiteRate.getText());
				HMMRateValues[Integer.parseInt(cmbxHMMCat.getSelectedItem().toString())-1]  = Double.parseDouble(txtHMMRate.getText());
				HMMProbValues[Integer.parseInt(cmbxHMMCat.getSelectedItem().toString())-1]  = Double.parseDouble(txtHMMProb.getText());

				inputvals = new DnaMLKData();
				inputvals.infile = txtInputFile.getText();
				inputvals.intree = txtInputTree.getText();
				inputvals.wgtsfile = txtWeightFile.getText();
				inputvals.catsfile = txtCatFile.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.outtree = txtOutTree.getText();
				inputvals.outtreeopt = "w";
				inputvals.TreeUseMethod = TreeSearchMethod.getSelectedItem().toString();
				inputvals.UseLengths = rdbtnUseLengthsYes.isSelected();
				inputvals.TTratio = Double.parseDouble(txtTTratio.getText());
				inputvals.useEmpBF = rdbtnEmpBFYes.isSelected();
				inputvals.BaseFreqA = Double.parseDouble(txtBaseFreqA.getText());
				inputvals.BaseFreqC = Double.parseDouble(txtBaseFreqC.getText());
				inputvals.BaseFreqG = Double.parseDouble(txtBaseFreqG.getText());
				inputvals.BaseFreqTU = Double.parseDouble(txtBaseFreqTU.getText());
				inputvals.OneCat = rdbtnCatYes.isSelected();
				inputvals.NumCats = Integer.parseInt(cmbxHMMcount.getSelectedItem().toString());
				
				inputvals.SiteRate1 = SiteRateValues[0];
				inputvals.SiteRate2 = SiteRateValues[1];
				inputvals.SiteRate3 = SiteRateValues[2];
				inputvals.SiteRate4 = SiteRateValues[3];
				inputvals.SiteRate5 = SiteRateValues[4];
				inputvals.SiteRate6 = SiteRateValues[5];
				inputvals.SiteRate7 = SiteRateValues[6];
				inputvals.SiteRate8 = SiteRateValues[7];
				inputvals.SiteRate9 = SiteRateValues[8];
				
				if (cmbxRateSite.getSelectedIndex() == 0)
				{
					inputvals.RateVar = "constant";
				}
				else if (cmbxRateSite.getSelectedIndex() == 1)
				{
					inputvals.RateVar = "gamma";
				}
				else if (cmbxRateSite.getSelectedIndex() == 2)
				{
					inputvals.RateVar = "invar";
				}
				else //(cmbxRateSite.getSelectedIndex() == 3)
				{
					inputvals.RateVar = "hmm";
				}
				//inputvals.RateVar = cmbxRateSite.getSelectedItem().toString();
				inputvals.AdjCor = rdbtnAdjCorYes.isSelected();
				inputvals.BlockLen = Double.parseDouble(txtBlockLen.getText());
				inputvals.CoeffVar = Double.parseDouble(txtCoeffVar.getText());
				inputvals.NumRates = Integer.parseInt(cmbxHMMcount.getSelectedItem().toString());
				
				inputvals.HMMrate1 = HMMRateValues[0];
				inputvals.HMMrate2 = HMMRateValues[1];
				inputvals.HMMrate3 = HMMRateValues[2];
				inputvals.HMMrate4 = HMMRateValues[3];
				inputvals.HMMrate5 = HMMRateValues[4];
				inputvals.HMMrate6 = HMMRateValues[5];
				inputvals.HMMrate7 = HMMRateValues[6];
				inputvals.HMMrate8 = HMMRateValues[7];
				inputvals.HMMrate9 = HMMRateValues[8];
				
				inputvals.HMMprob1 = HMMProbValues[0];
				inputvals.HMMprob2 = HMMProbValues[1];
				inputvals.HMMprob3 = HMMProbValues[2];
				inputvals.HMMprob4 = HMMProbValues[3];
				inputvals.HMMprob5 = HMMProbValues[4];
				inputvals.HMMprob6 = HMMProbValues[5];
				inputvals.HMMprob7 = HMMProbValues[6];
				inputvals.HMMprob8 = HMMProbValues[7];
				inputvals.HMMprob9 = HMMProbValues[8];
				
				inputvals.InvarFract = Double.parseDouble(txtFracInvar.getText());
				inputvals.SitesWeight = rdbtnSitesWeightYes.isSelected();
				inputvals.GlobalRe = rdbtnGlobalReYes.isSelected();
				inputvals.RandInput = rdbtnRandOrderYes.isSelected();
				inputvals.RandNum = Integer.parseInt(txtRandSeed.getText());
				inputvals.Njumble = Integer.parseInt(txtNumberJumble.getText());
				inputvals.MultData = rdbtnAnalyzeMultDataYes.isSelected();
				inputvals.MultDSet = rdbtnDataSets.isSelected();
				inputvals.NumSeqs = Integer.parseInt(txtNumSeqs.getText());
				inputvals.InputSeq = rdbtnInputSeqYes.isSelected();
				inputvals.PrintData = rdbtnPrintDataYes.isSelected();
				inputvals.PrintInd = rdbtnPrintIndYes.isSelected();
				inputvals.PrintTree = rdbtnPrintTreeYes.isSelected();
				inputvals.WriteTree = rdbtnWriteTreeYes.isSelected();
				inputvals.DotDiff = rdbtnDotDiffYes.isSelected();
				inputvals.RecHypo = rdbtnRecHypoYes.isSelected();

				btnExecute.setEnabled(false);	
				String title = "Dnamlk Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread dnaMLKThread = new Thread() {
						public void run() {
							runDnaMLKThreads();
						}
			  	    };
			  	    dnaMLKThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 4 21,alignx right");
		
		btnQuit = new JButton("Quit");
		btnQuit.setBounds(975, 576, 89, 23);
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				if(phylipCall)
				{
					frmDnaMLKControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 5 21");
	}
	
	public boolean checkInputVals(){
		TestFileNames test = new TestFileNames();
		
		if (!test.DuplicateFileNames(inputvals.infile, "Input", inputvals.outfile, "Output"))
		{			
			return false;		
		}
		
		if (!test.DuplicateFileNames(inputvals.intree, "Intree", inputvals.outtree, "Outtree"))
		{			
			return false;		
		}

		if (!test.FileAvailable(inputvals.infile, "Input"))
		{
			return false;
		}
		
		if (inputvals.SitesWeight){
			if (!test.FileAvailable(inputvals.wgtsfile, "Weights"))
			{
				return false;
			}
		}
		
		if (!inputvals.OneCat){
			if (!test.FileAvailable(inputvals.catsfile, "Categories"))
			{
				return false;
			}
		}

		String opt = test.FileAlreadyExists(inputvals.outfile, "Outfile");
		if (opt == "q")
		{
			return false;
		}
		else
		{
			inputvals.outfileopt = opt;
		}

		opt = test.FileAlreadyExists(inputvals.outtree, "Outtree");
		if (opt == "q")
		{
			return false;
		}
		else
		{
			inputvals.outfileopt = opt;
		}


		if (inputvals.BaseFreqA < 0) {
			String msg1 = "Input value: Base Frequency for A cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.BaseFreqC < 0) {
			String msg1 = "Input value: Base Frequency for C cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.BaseFreqG < 0) {
			String msg1 = "Input value: Base Frequency for G cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.BaseFreqTU < 0) {
			String msg1 = "Input value: Base Frequency for T/U cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		double sum = inputvals.BaseFreqA + inputvals.BaseFreqC + inputvals.BaseFreqG + inputvals.BaseFreqTU;
		inputvals.BaseFreqA = inputvals.BaseFreqA / sum;
		inputvals.BaseFreqC = inputvals.BaseFreqC / sum;
		inputvals.BaseFreqG = inputvals.BaseFreqG / sum;
		inputvals.BaseFreqTU = inputvals.BaseFreqTU / sum;
		
		if (inputvals.BlockLen <= 0) {
			String msg1 = "Input value: Mean Block Length must be greater than 0.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.TTratio < 0) {
			String msg1 = "Input value: Transition / Transversion Ratio cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.CoeffVar <= 0.0) {
			String msg1 = "Input value: Coefficient of Variation must be greater than 0.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if ((inputvals.RandNum % 2) == 0)
		{
			String msg1 = "Random number seed must be odd.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;			
		}

		// autoscale the probabilities 
		// unrolled because of JNA limitations (BLEH!)
		int i = 1;
		double probsum = inputvals.HMMprob1;
		
		if (i<inputvals.NumRates){
   		probsum += inputvals.HMMprob2;
		}
		else
		{
			inputvals.HMMprob2 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumRates){
   		probsum += inputvals.HMMprob3;
		}
		else
		{
			inputvals.HMMprob3 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumRates){
   		probsum += inputvals.HMMprob4;
		}
		else
		{
			inputvals.HMMprob4 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumRates){
   		probsum += inputvals.HMMprob5;
		}
		else
		{
			inputvals.HMMprob5 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumRates){
   		probsum += inputvals.HMMprob6;
		}
		else
		{
			inputvals.HMMprob6 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumRates){
   		probsum += inputvals.HMMprob7;
		}
		else
		{
			inputvals.HMMprob7 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumRates){
   		probsum += inputvals.HMMprob8;
		}
		else
		{
			inputvals.HMMprob8 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumRates){
   		probsum += inputvals.HMMprob9;
		}
		else
		{
			inputvals.HMMprob9 = 0.0;
		}
		
		inputvals.HMMprob1 = inputvals.HMMprob1/probsum;
		inputvals.HMMprob2 = inputvals.HMMprob2/probsum;
		inputvals.HMMprob3 = inputvals.HMMprob3/probsum;
		inputvals.HMMprob4 = inputvals.HMMprob4/probsum;
		inputvals.HMMprob5 = inputvals.HMMprob5/probsum;
		inputvals.HMMprob6 = inputvals.HMMprob6/probsum;
		inputvals.HMMprob7 = inputvals.HMMprob7/probsum;
		inputvals.HMMprob8 = inputvals.HMMprob8/probsum;
		inputvals.HMMprob9 = inputvals.HMMprob9/probsum;
		
		return true;
	}
	
	protected void runDnaMLKThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("dnamlk", DnaMLK.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("DnaMLK");
    		return;
    	}
		try 
		{
	  	    Thread dnaMLKRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			DnaMLK dnamlk = (DnaMLK) Native.loadLibrary("dnamlk", DnaMLK.class);
		  	        dnamlk.dnamlk(		
		  	       		inputvals.infile,
		  	       		inputvals.intree,
		  	       		inputvals.wgtsfile,
		  	       		inputvals.catsfile,
		  	       		inputvals.outfile,
		  	       		inputvals.outfileopt,
		  	       		inputvals.outtree,
		  	       		inputvals.outtreeopt,
		  	       		inputvals.TreeUseMethod,
		  	       		inputvals.UseLengths,
		  	       		inputvals.TTratio,
		  	       		inputvals.useEmpBF,
		  	       		inputvals.BaseFreqA,
		  	       		inputvals.BaseFreqC,
		  	       		inputvals.BaseFreqG,
		  	       		inputvals.BaseFreqTU,
		  	       		inputvals.OneCat,
		  	       		inputvals.NumCats,
		  	       		inputvals.SiteRate1,
		  				inputvals.SiteRate2,
		  				inputvals.SiteRate3,
		  				inputvals.SiteRate4,
		  				inputvals.SiteRate5,
		  				inputvals.SiteRate6,
		  				inputvals.SiteRate7,
		  				inputvals.SiteRate8,
		  				inputvals.SiteRate9,
		  	       		inputvals.RateVar,
		  	       		inputvals.AdjCor,
		  	       		inputvals.BlockLen,
		  	       		inputvals.CoeffVar,
		  	       		inputvals.NumRates,
		  	       		inputvals.HMMrate1,
		  	            inputvals.HMMrate2,
		  	            inputvals.HMMrate3,
		  	            inputvals.HMMrate4,
		  	            inputvals.HMMrate5,
		  	            inputvals.HMMrate6,
		  	            inputvals.HMMrate7,
		  	            inputvals.HMMrate8,
		  	            inputvals.HMMrate9,
		  	       		inputvals.HMMprob1,
		  	            inputvals.HMMprob2,
		  	            inputvals.HMMprob3,
		  	            inputvals.HMMprob4,
		  	            inputvals.HMMprob5,
		  	            inputvals.HMMprob6,
		  	            inputvals.HMMprob7,
		  	            inputvals.HMMprob8,
		  	            inputvals.HMMprob9,
		  	            inputvals.InvarFract,
		  	           	inputvals.SitesWeight,
		  	       		inputvals.GlobalRe,
		  	       		inputvals.RandInput,
		  	       		inputvals.RandNum,
		  	       		inputvals.Njumble,
		  	       		inputvals.MultData,
		  	       		inputvals.MultDSet,
		  	       		inputvals.NumSeqs,
		  	       		inputvals.InputSeq,
		  	       		inputvals.PrintData,
		  	       		inputvals.PrintInd,
		  	       		inputvals.PrintTree,
		  	       		inputvals.WriteTree,
		  	       		inputvals.DotDiff,
		  	       		inputvals.RecHypo); 
		  	    };
	  	    };
	  	  dnaMLKRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (dnaMLKRunThread.isAlive());
	  	    }
		} 
		catch (InterruptedException e) {
			if (inputvals.PrintInd)
			{
				updateProgress();
			}
		}
	}
	  
	private void updateProgress(){
		SwingUtilities.invokeLater(new Runnable(){
			public void run() {
				DisplayProgress oldp = dp;
				dp = new DisplayProgress(inTitle,inCurdir);
				// doing the dispose in this order keeps the main window from flickering
				if (oldp != null) 
				{
					oldp.dispose();
				}
			}
		});
	}
}
