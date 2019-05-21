package phylip;

import java.awt.EventQueue;
import javax.swing.JFileChooser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

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
import utilities.CopyFile;

import com.sun.jna.Library;
import com.sun.jna.Native;

import consense.ConsenseUserInterface;

import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

import drawgram.DrawgramUserInterface;
import drawtree.DrawtreeUserInterface;


public class DnaMLUserInterface {
   public interface DnaML extends Library {
        public void dnaml(
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
        		boolean SpeedAn,
        		boolean GlobalRe,
        		boolean RandInput,
        		int RandNum,
        		int Njumble,
        		boolean OutRoot,
        		int OutNum,
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

	public class DnaMLData {
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
		boolean SpeedAn;
		boolean GlobalRe;
		boolean RandInput;
		int RandNum;
		int Njumble;
		boolean OutRoot;
		int OutNum;
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

	private DnaMLData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean ExplicitWgts;
	private boolean phylipCall;
	private boolean bootstrapCall;
	private String initFile;

	private double[] SiteRateValues;
	private double[] HMMRateValues;
	private double[] HMMProbValues;
	private int lastsiteratecat;
	private int lasthmmcat;

	private JFrame frmDnaMLControls;
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
	private JLabel lblSpeedAn;
	private JRadioButton rdbtnSpeedAnYes;
	private JRadioButton rdbtnSpeedAnNo;
	private JLabel lblGlobalRe;
	private JRadioButton rdbtnGlobalReYes;
	private JRadioButton rdbtnGlobalReNo;
	private JLabel lblRandOrder;
	private JRadioButton rdbtnRandOrderYes;
	private JRadioButton rdbtnRandOrderNo;
	private JLabel lblRandSeed;
	private JTextField txtRandSeed;
	private JLabel lblOutRoot;
	private JRadioButton rdbtnOutRootYes;
	private JRadioButton rdbtnOutRootNo;
	private JLabel lblOutRootNum;
	private JTextField txtOutRootNum;
	private JLabel lblAnalyzeMultData;
	private JRadioButton rdbtnAnalyzeMultDataYes;
	private JRadioButton rdbtnAnalyzeMultDataNo;
	private JLabel lblMultData;
	private JRadioButton rdbtnDataSets;
	private JRadioButton rdbtnWeights;
	private JLabel lblInputSeq;
	private JRadioButton rdbtnInputFileSeqYes;
	private JRadioButton rdbtnInputFileSeqNo;
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
	private JComboBox cmbxRateCatnum;
	private JComboBox cmbxTreeSearchMethod;
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
	private JButton btnDefaults;
	private JButton btnStored;
	
	private JComboBox cmbxDisplayKind;
	private JButton btnDisplayTree;


	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DnaMLUserInterface window = new DnaMLUserInterface(args);
					window.frmDnaMLControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmDnaMLControls.getRootPane());
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
			cmbxRateCatnum.setEnabled(false);
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
			cmbxRateCatnum.setEnabled(true);
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
	
	protected void SpeedAnToggle(boolean isSpeedAn){
		if(isSpeedAn){
			rdbtnSpeedAnYes.setSelected(true);
			rdbtnSpeedAnNo.setSelected(false);
		}
		else{
			rdbtnSpeedAnYes.setSelected(false);
			rdbtnSpeedAnNo.setSelected(true);
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
	
	protected void OutRootToggle(boolean isOut){
		if(isOut){
			rdbtnOutRootNo.setSelected(true);
			rdbtnOutRootYes.setSelected(false);
			lblOutRootNum.setEnabled(false);
			txtOutRootNum.setEnabled(false);
		}
		else{
			rdbtnOutRootNo.setSelected(false);
			rdbtnOutRootYes.setSelected(true);
			lblOutRootNum.setEnabled(true);
			txtOutRootNum.setEnabled(true);
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
			rdbtnInputFileSeqYes.setEnabled(true);
			rdbtnInputFileSeqNo.setEnabled(true);
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
			rdbtnInputFileSeqYes.setEnabled(false);
			rdbtnInputFileSeqNo.setEnabled(false);
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
			rdbtnInputFileSeqYes.setSelected(true);
			rdbtnInputFileSeqNo.setSelected(false);
		}
		else{
			rdbtnInputFileSeqYes.setSelected(false);
			rdbtnInputFileSeqNo.setSelected(true);
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
			lblRandOrder.setEnabled(true);
			rdbtnRandOrderYes.setEnabled(true);
			rdbtnRandOrderNo.setEnabled(true);
			RandOrderToggle(rdbtnRandOrderYes.isSelected());
		}
		else{
			txtInputTree.setEnabled(true);
			btnInputTree.setEnabled(true);
			lblUseLengths.setEnabled(true);
			rdbtnUseLengthsYes.setEnabled(true);
			rdbtnUseLengthsNo.setEnabled(true);
			lblRandOrder.setEnabled(false);
			rdbtnRandOrderYes.setEnabled(false);
			rdbtnRandOrderNo.setEnabled(false);
			lblRandSeed.setEnabled(false);
			txtRandSeed.setEnabled(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEnabled(false);
			lblRandOdd.setEnabled(false);
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
	
	protected void CallPlotType(int plotkind){
		
		if (txtOutTree.getText() != "DnaMLOuttree")
		{
			// copy outtree to DnaMLOuttree if it isn't already there
			try {
				new CopyFile(txtOutTree.getText(), "DnaMLOuttree");
	        } catch ( IOException e ) {
	             e.printStackTrace();
	        }
		}
		
		final String[] plottree = new String[] {"DnaML"};
		if(plotkind == 0) // Drawgram
		{
			try {
				DrawgramUserInterface.main(plottree);
			} catch (Exception err) {
				err.printStackTrace();
			}						
		}
		else //if (plotkind == 1) // Drawtree
		{
			try {
				DrawtreeUserInterface.main(plottree);
			} catch (Exception err) {
				err.printStackTrace();
			}						
		}
	}
	
	protected void DoInit()
	{
		// reset everything if there is an init file
		getStoredSettings();
	}

	/**
	 * Create the application.
	 */
	public DnaMLUserInterface(String[] args) {
		phylipCall = false;
		bootstrapCall = false;
		if (args.length > 0)
		{
			if (args[0].contains("Phylip"))
			{
				phylipCall = true;
			}
			if (args.length > 1)
			{
				if (args[1].contains("bootstrap"))
				{
					bootstrapCall = true;
				}
			}
		}
		initialize();
		DoInit();
		
		if (bootstrapCall){
  			txtOutTree.setText("DnaMLOuttree");
  			getSeqBootSettings();
			MultToggle(true);
			DataWeightToggle(true);
 			}
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
		
		initFile = "dnamlInit.txt";

		frmDnaMLControls = new JFrame();
		frmDnaMLControls.setBackground(new Color(204, 255, 255));
		frmDnaMLControls.setTitle("Dnaml");
		frmDnaMLControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDnaMLControls.setBounds(100, 100, 1200, 810);
		frmDnaMLControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDnaMLControls.setPreferredSize(new Dimension(frmDnaMLControls.getBounds().width, frmDnaMLControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmDnaMLControls.getPreferredSize());
		frmDnaMLControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmDnaMLControls.getPreferredSize());
		
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow][pref!,grow][pref!,grow][pref!,grow]", "[][][][][][][][][][][][][][][][][][][][][][][][]"));
	
		btnInputFile = new JButton("Input File");
		btnInputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputFile);
			}
		});
		panel.add(btnInputFile, "cell 0 0,growx");
		
		txtInputFile = new JTextField();
		txtInputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInputFile.setBounds(168, 11, 925, 20);
		panel.add(txtInputFile, "cell 1 0 5 1,growx");
		
		btnInputTree = new JButton("Input Tree");
		btnInputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputTree);
			}
		});
		btnInputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnInputTree, "cell 0 1,growx");
		
		txtInputTree = new JTextField();
		txtInputTree.setEnabled(false);
		txtInputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInputTree, "cell 1 1 5 1,growx");
	
		btnWeightFile = new JButton("Weights File");
		btnWeightFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightFile);
			}
		});
		btnWeightFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnWeightFile, "cell 0 2,growx");
	
		txtWeightFile = new JTextField();
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtWeightFile, "cell 1 2 5 1,growx");
		
		btnCatFile = new JButton("Categories File");
		btnCatFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtCatFile);
			}
		});
		btnCatFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnCatFile, "cell 0 3,growx");
		
		txtCatFile = new JTextField();
		txtCatFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtCatFile, "cell 1 3 5 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		panel.add(btnOutputFile, "cell 0 4,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutputFile, "cell 1 4 5 1,growx");

		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutTree);
			}
		});
		panel.add(btnOutputTree, "cell 0 5,growx");
		
		txtOutTree = new JTextField();
		txtOutTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutTree, "cell 1 5 5 1,growx");
		
		lblSearchBest = new JLabel("Search for best tree:");
		lblSearchBest.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSearchBest.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSearchBest, "flowx,cell 0 6 2 1,alignx right");
		
		cmbxTreeSearchMethod = new JComboBox();
		cmbxTreeSearchMethod.setModel(new DefaultComboBoxModel(new String[] {"Yes", "No, use user trees in input file", "Yes, rearrange on user tree"}));
		cmbxTreeSearchMethod.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(cmbxTreeSearchMethod.getSelectedIndex());
			}
		});
		panel.add(cmbxTreeSearchMethod, "cell 2 6 2 1,growx");
		
		lblUseLengths = new JLabel("Use lengths from user trees:");
		lblUseLengths.setHorizontalAlignment(SwingConstants.RIGHT);
		lblUseLengths.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseLengths, "flowx,cell 0 7 2 1,alignx right");
		
		rdbtnUseLengthsYes = new JRadioButton("Yes");
		rdbtnUseLengthsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseLengthsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UseLengthsToggle(true);
			}
		});
		rdbtnUseLengthsYes.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnUseLengthsYes, "cell 2 7");
		
		rdbtnUseLengthsNo = new JRadioButton("No");
		rdbtnUseLengthsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseLengthsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UseLengthsToggle(false);
			}
		});
		rdbtnUseLengthsNo.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnUseLengthsNo, "cell 2 7");

		lblTTratio = new JLabel("Transition/transversion ratio:");
		lblTTratio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblTTratio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTTratio, "flowx,cell 0 8 2 1,alignx right");
		
		txtTTratio = new JTextField();
		txtTTratio.setFont(new Font("Arial", Font.PLAIN, 13));
		txtTTratio.setColumns(6);
		panel.add(txtTTratio, "cell 2 8");
		
		lblEmpBF = new JLabel("Use empirical base frequencies:");
		lblEmpBF.setHorizontalAlignment(SwingConstants.RIGHT);
		lblEmpBF.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblEmpBF, "flowx,cell 0 9 2 1,alignx right");
		
		rdbtnEmpBFYes = new JRadioButton("Yes");
		rdbtnEmpBFYes.setBackground(new Color(204, 255, 255));
		rdbtnEmpBFYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnEmpBFYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EmpBFToggle(true);
			}
		});
		panel.add(rdbtnEmpBFYes, "cell 2 9");
		
		rdbtnEmpBFNo = new JRadioButton("No");
		rdbtnEmpBFNo.setBackground(new Color(204, 255, 255));
		rdbtnEmpBFNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnEmpBFNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EmpBFToggle(false);
			}
		});
		panel.add(rdbtnEmpBFNo, "cell 2 9");
		
		lblBaseFreq = new JLabel("Base frequencies:");
		lblBaseFreq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblBaseFreq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreq, "flowx,cell 0 10 2 1,alignx right");
		
		lblBaseFreqA = new JLabel("A:");
		lblBaseFreqA.setBackground(new Color(153, 255, 255));
		lblBaseFreqA.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqA, "cell 2 10");
		
		txtBaseFreqA = new JTextField();
		txtBaseFreqA.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqA.setColumns(6);
		panel.add(txtBaseFreqA, "cell 2 10");
		
		lblBaseFreqC = new JLabel("   C:");
		lblBaseFreqC.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqC, "cell 3 10");
		
		txtBaseFreqC = new JTextField();
		txtBaseFreqC.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqC.setColumns(6);
		panel.add(txtBaseFreqC, "cell 3 10");
		
		lblBaseFreqG = new JLabel("G:");
		lblBaseFreqG.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqG, "cell 2 11");
		
		txtBaseFreqG = new JTextField();
		txtBaseFreqG.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqG.setColumns(6);
		panel.add(txtBaseFreqG, "cell 2 11");
		
		lblBaseFreqTU = new JLabel("T/U:");
		lblBaseFreqTU.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqTU, "cell 3 11");
		
		txtBaseFreqTU = new JTextField();
		txtBaseFreqTU.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqTU.setColumns(6);
		panel.add(txtBaseFreqTU, "cell 3 11");
		
		lblCatSites = new JLabel("One category of sites:");
		lblCatSites.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCatSites.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCatSites, "flowx,cell 0 12 2 1,alignx right");
		
		rdbtnCatYes = new JRadioButton("Yes");
		rdbtnCatYes.setBackground(new Color(204, 255, 255));
		rdbtnCatYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCatYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CatToggle(true);
			}
		});
		panel.add(rdbtnCatYes, "cell 2 12");
		
		rdbtnCatNo = new JRadioButton("No");
		rdbtnCatNo.setBackground(new Color(204, 255, 255));
		rdbtnCatNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCatNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CatToggle(false);
			}
		});
		panel.add(rdbtnCatNo, "cell 2 12");
		
		lblNumCat = new JLabel("Number of site categories:");
		lblNumCat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumCat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNumCat, "flowx,cell 0 13 2 1,alignx right");
		
		cmbxNumCat = new JComboBox();
		cmbxNumCat.setEnabled(false);
		cmbxNumCat.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxNumCat.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxNumCat, "cell 2 13");
		
		lblcatnum = new JLabel("Category:");
		lblcatnum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblcatnum.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblcatnum, "cell 2 13");
		
		cmbxRateCatnum = new JComboBox();
		cmbxRateCatnum.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxRateCatnum.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxRateCatnum.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplaySiteRateValue(Integer.parseInt(cmbxRateCatnum.getSelectedItem().toString()));
			}
		});
		panel.add(cmbxRateCatnum, "cell 3 13");
		
		lblRate = new JLabel("Rate:");
		lblRate.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRate.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRate, "cell 3 13");
		
		txtSiteRate = new JTextField();
		txtSiteRate.setHorizontalAlignment(SwingConstants.CENTER);
		txtSiteRate.setFont(new Font("Arial", Font.PLAIN, 13));
		txtSiteRate.setColumns(5);
		panel.add(txtSiteRate, "cell 3 13");
		
		lblRateSite = new JLabel("Rate variation among sites:");
		lblRateSite.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRateSite.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRateSite, "flowx,cell 0 14 2 1,alignx right");
		
		cmbxRateSite = new JComboBox();
		cmbxRateSite.setModel(new DefaultComboBoxModel(new String[] {"Constant rate", "Gamma distance rates", "Gamma + invariant sites","User-defined HMM of rates"}));
		cmbxRateSite.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxRateSite.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				RatesiteToggle(cmbxRateSite.getSelectedIndex());
			}
		});
		panel.add(cmbxRateSite, "cell 2 14 2 1,growx");
		
		lblRatesAdjCor = new JLabel("Rates on adjacent sites correlated:");
		lblRatesAdjCor.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRatesAdjCor.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRatesAdjCor, "flowx,cell 0 15 2 1,alignx right");
		
		rdbtnAdjCorYes = new JRadioButton("Yes");
		rdbtnAdjCorYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAdjCorYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AdjCorToggle(true);
			}
		});
		rdbtnAdjCorYes.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnAdjCorYes, "cell 2 15");
		
		rdbtnAdjCorNo = new JRadioButton("No");
		rdbtnAdjCorNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAdjCorNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AdjCorToggle(false);
			}
		});
		rdbtnAdjCorNo.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnAdjCorNo, "cell 2 15");
		
		lblMeanLen = new JLabel("Mean block length:");
		lblMeanLen.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMeanLen.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMeanLen, "cell 3 15");
		
		txtBlockLen = new JTextField();
		txtBlockLen.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBlockLen.setColumns(5);
		panel.add(txtBlockLen, "cell 3 15");
		
		lblCoeffVar = new JLabel("Coefficient of variation:");
		lblCoeffVar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCoeffVar.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCoeffVar, "flowx,cell 0 16 2 1,alignx right");

		txtCoeffVar = new JTextField();
		txtCoeffVar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCoeffVar.setColumns(6);
		panel.add(txtCoeffVar, "cell 2 16");
		
		lblCoeffVarNote = new JLabel("(for gamma dist = 1/\u221Aalpha)");
		lblCoeffVarNote.setHorizontalAlignment(SwingConstants.LEFT);
		lblCoeffVarNote.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCoeffVarNote, "cell 3 16");
		
		lblHMM = new JLabel("Rates for HMM:");
		lblHMM.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMM.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblHMM, "flowx,cell 0 17 2 1,alignx right");
		
		lblInvar = new JLabel("");
		lblInvar.setHorizontalAlignment(SwingConstants.LEFT);
		lblInvar.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblInvar, "cell 2 17");
		
		lblNumHMMcats = new JLabel("Number of rate categories:");
		lblNumHMMcats.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumHMMcats.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNumHMMcats, "flowx,cell 0 18 2 1,alignx right");

		cmbxHMMcount = new JComboBox();
		cmbxHMMcount.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxHMMcount.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxHMMcount, "cell 2 18");
		
		lblHMMCat = new JLabel("Category:");
		lblHMMCat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMMCat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblHMMCat, "cell 3 18");
		
		cmbxHMMCat = new JComboBox();
		cmbxHMMCat.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxHMMCat.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxHMMCat.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplayHMMValues(Integer.parseInt(cmbxHMMCat.getSelectedItem().toString()));
			}
		});
		panel.add(cmbxHMMCat, "cell 3 18");
		
		lblHMMRate = new JLabel("Rate:");
		lblHMMRate.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMMRate.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblHMMRate, "cell 2 19");
			
		txtHMMRate = new JTextField();
		txtHMMRate.setHorizontalAlignment(SwingConstants.CENTER);
		txtHMMRate.setFont(new Font("Arial", Font.PLAIN, 13));
		txtHMMRate.setColumns(5);
		panel.add(txtHMMRate, "cell 2 19");
		
		lblProb = new JLabel("Probablilty:");
		lblProb.setHorizontalAlignment(SwingConstants.RIGHT);
		lblProb.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblProb, "cell 3 19");
		
		txtHMMProb = new JTextField();
		txtHMMProb.setHorizontalAlignment(SwingConstants.CENTER);
		txtHMMProb.setFont(new Font("Arial", Font.PLAIN, 13));
		txtHMMProb.setColumns(5);
		panel.add(txtHMMProb, "cell 3 19");
		
		lblFracInvar = new JLabel("Fraction of sites invariant:");
		lblFracInvar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblFracInvar.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblFracInvar, "flowx,cell 0 20 2 1,alignx right");
	
		txtFracInvar = new JTextField();
		txtFracInvar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtFracInvar.setColumns(10);
		panel.add(txtFracInvar, "cell 2 20");

		lblSitesWeight = new JLabel("Sites weighted:");
		lblSitesWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSitesWeight.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSitesWeight, "flowx,cell 0 21 2 1,alignx right");
		
		rdbtnSitesWeightYes = new JRadioButton("Yes");
		rdbtnSitesWeightYes.setBackground(new Color(204, 255, 255));
		rdbtnSitesWeightYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesWeightYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SitesWeightToggle(false);
			}
		});
		panel.add(rdbtnSitesWeightYes, "cell 2 21");
		
		rdbtnSitesWeightNo = new JRadioButton("No");
		rdbtnSitesWeightNo.setBackground(new Color(204, 255, 255));
		rdbtnSitesWeightNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesWeightNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SitesWeightToggle(true);
			}
		});
		panel.add(rdbtnSitesWeightNo, "cell 2 21");
		
		lblSpeedAn = new JLabel("Speedier but rougher analysis:");
		lblSpeedAn.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSpeedAn.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSpeedAn, "flowx,cell 0 22 2 1,alignx right");
		
		rdbtnSpeedAnYes = new JRadioButton("Yes");
		rdbtnSpeedAnYes.setBackground(new Color(204, 255, 255));
		rdbtnSpeedAnYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSpeedAnYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SpeedAnToggle(true);
			}
		});
		panel.add(rdbtnSpeedAnYes, "cell 2 22");
		
		rdbtnSpeedAnNo = new JRadioButton("No");
		rdbtnSpeedAnNo.setBackground(new Color(204, 255, 255));
		rdbtnSpeedAnNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSpeedAnNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SpeedAnToggle(false);
			}
		});
		panel.add(rdbtnSpeedAnNo, "cell 2 22");
		
		// **** column 2 start ****
		lblGlobalRe = new JLabel("Global rearrangements:");
		lblGlobalRe.setHorizontalAlignment(SwingConstants.RIGHT);
		lblGlobalRe.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGlobalRe, "cell 4 6,alignx right");
		
		rdbtnGlobalReYes = new JRadioButton("Yes");
		rdbtnGlobalReYes.setBackground(new Color(204, 255, 255));
		rdbtnGlobalReYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnGlobalReYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalReToggle(false);
			}
		});
		panel.add(rdbtnGlobalReYes, "cell 5 6");
		
		rdbtnGlobalReNo = new JRadioButton("No");
		rdbtnGlobalReNo.setBackground(new Color(204, 255, 255));
		rdbtnGlobalReNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnGlobalReNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalReToggle(true);
			}
		});
		panel.add(rdbtnGlobalReNo, "cell 5 6");
		
		lblRandOrder = new JLabel("Randomize input order of sequences:");
		lblRandOrder.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandOrder.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOrder, "cell 4 7,alignx right");
		
		rdbtnRandOrderYes = new JRadioButton("Yes");
		rdbtnRandOrderYes.setBackground(new Color(204, 255, 255));
		rdbtnRandOrderYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRandOrderYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(true);
			}
		});
		panel.add(rdbtnRandOrderYes, "cell 5 7");
		
		rdbtnRandOrderNo = new JRadioButton("No, use input order");
		rdbtnRandOrderNo.setBackground(new Color(204, 255, 255));
		rdbtnRandOrderNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRandOrderNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(false);
			}
		});
		panel.add(rdbtnRandOrderNo, "cell 5 7");
		
		lblRandSeed = new JLabel("Random number seed:");
		lblRandSeed.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandSeed.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandSeed, "cell 4 8,alignx right");
		
		txtRandSeed = new JTextField();
		txtRandSeed.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRandSeed.setColumns(6);
		panel.add(txtRandSeed, "cell 5 8");
		
		lblRandOdd = new JLabel("(must be odd)");
		lblRandOdd.setHorizontalAlignment(SwingConstants.LEFT);
		lblRandOdd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOdd, "cell 5 8");

		lblNumberJumble = new JLabel("Number of times to jumble:");
		lblNumberJumble.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberJumble.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblNumberJumble, "cell 4 9,alignx right");

		txtNumberJumble = new JTextField();
		txtNumberJumble.setColumns(6);
		txtNumberJumble.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberJumble.setBounds(826, 248, 86, 20);
		panel.add(txtNumberJumble, "cell 5 9");
	
		lblOutRoot = new JLabel("Outgroup root:");
		lblOutRoot.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutRoot.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblOutRoot, "cell 4 10,alignx right");
		
		rdbtnOutRootYes = new JRadioButton("Yes");
		rdbtnOutRootYes.setBackground(new Color(204, 255, 255));
		rdbtnOutRootYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutRootYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutRootToggle(false);
			}
		});
		panel.add(rdbtnOutRootYes, "cell 5 10");
		
		rdbtnOutRootNo = new JRadioButton("No, use as outgroup species");
		rdbtnOutRootNo.setBackground(new Color(204, 255, 255));
		rdbtnOutRootNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutRootNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutRootToggle(true);
			}
		});
		panel.add(rdbtnOutRootNo, "cell 5 10");
		
		lblOutRootNum = new JLabel("Number of the outgroup:");
		lblOutRootNum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutRootNum.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblOutRootNum, "cell 4 11,alignx right");
		
		txtOutRootNum = new JTextField();
		txtOutRootNum.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutRootNum.setColumns(6);
		panel.add(txtOutRootNum, "cell 5 11");
		
		lblAnalyzeMultData = new JLabel("Analyze multiple data sets:");
		lblAnalyzeMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAnalyzeMultData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAnalyzeMultData, "cell 4 12,alignx right");
		
		rdbtnAnalyzeMultDataYes = new JRadioButton("Yes");
		rdbtnAnalyzeMultDataYes.setBackground(new Color(204, 255, 255));
		rdbtnAnalyzeMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAnalyzeMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		panel.add(rdbtnAnalyzeMultDataYes, "cell 5 12");
		
		rdbtnAnalyzeMultDataNo = new JRadioButton("No");
		rdbtnAnalyzeMultDataNo.setBackground(new Color(204, 255, 255));
		rdbtnAnalyzeMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAnalyzeMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		panel.add(rdbtnAnalyzeMultDataNo, "cell 5 12");
		
		lblMultData = new JLabel("Multiple data sets or multiple weights:");
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMultData, "cell 4 13,alignx right");
		
		rdbtnDataSets = new JRadioButton("Data sets");
		rdbtnDataSets.setBackground(new Color(204, 255, 255));
		rdbtnDataSets.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(true);
			}
		});
		rdbtnDataSets.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(rdbtnDataSets, "cell 5 13");
		
		rdbtnWeights = new JRadioButton("Weights");
		rdbtnWeights.setBackground(new Color(204, 255, 255));
		rdbtnWeights.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(false);
			}
		});		
		rdbtnWeights.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(rdbtnWeights, "cell 5 13");

		lblHowManyData = new JLabel("Number:");
		lblHowManyData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHowManyData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblHowManyData, "cell 5 14");
		
		txtNumSeqs = new JTextField();
		txtNumSeqs.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumSeqs.setColumns(6);
		panel.add(txtNumSeqs, "cell 5 14");

		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblInputSeq, "cell 4 15,alignx right");
		
		rdbtnInputFileSeqYes = new JRadioButton("Interleaved");
		rdbtnInputFileSeqYes.setBackground(new Color(204, 255, 255));
		rdbtnInputFileSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputFileSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		panel.add(rdbtnInputFileSeqYes, "cell 5 15");
		
		rdbtnInputFileSeqNo = new JRadioButton("Sequential");
		rdbtnInputFileSeqNo.setBackground(new Color(204, 255, 255));
		rdbtnInputFileSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputFileSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		panel.add(rdbtnInputFileSeqNo, "cell 5 15");
		
		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintData, "cell 4 16,alignx right");
		
		rdbtnPrintDataYes = new JRadioButton("Yes");
		rdbtnPrintDataYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		panel.add(rdbtnPrintDataYes, "cell 5 16");
		
		rdbtnPrintDataNo = new JRadioButton("No");
		rdbtnPrintDataNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		panel.add(rdbtnPrintDataNo, "cell 5 16");
		
		lblPrintTree = new JLabel("Print out tree:");
		lblPrintTree.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintTree, "cell 4 17,alignx right");
		
		rdbtnPrintTreeYes = new JRadioButton("Yes");
		rdbtnPrintTreeYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(true);
			}
		});
		panel.add(rdbtnPrintTreeYes, "cell 5 17");
		
		rdbtnPrintTreeNo = new JRadioButton("No");
		rdbtnPrintTreeNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(false);
			}
		});
		panel.add(rdbtnPrintTreeNo, "cell 5 17");
		
		lblWriteTree = new JLabel("Write out trees onto tree file:");
		lblWriteTree.setHorizontalAlignment(SwingConstants.RIGHT);
		lblWriteTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblWriteTree, "cell 4 18,alignx right");
		
		rdbtnWriteTreeYes = new JRadioButton("Yes");
		rdbtnWriteTreeYes.setBackground(new Color(204, 255, 255));
		rdbtnWriteTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(true);
			}
		});
		panel.add(rdbtnWriteTreeYes, "cell 5 18");
		
		rdbtnWriteTreeNo = new JRadioButton("No");
		rdbtnWriteTreeNo.setBackground(new Color(204, 255, 255));
		rdbtnWriteTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(false);
			}
		});
		panel.add(rdbtnWriteTreeNo, "cell 5 18");

		lblDotDiff = new JLabel("Use dot-differencing to display them:");
		lblDotDiff.setFont(new Font("Arial", Font.BOLD, 13));
		lblDotDiff.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblDotDiff, "cell 4 19,alignx right");

		rdbtnDotDiffYes = new JRadioButton("Yes");
		rdbtnDotDiffYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(true);
			}
		});
		panel.add(rdbtnDotDiffYes, "cell 5 19");

		rdbtnDotDiffNo = new JRadioButton("No");
		rdbtnDotDiffNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(false);
			}
		});
		panel.add(rdbtnDotDiffNo, "cell 5 19");
		
		lblRecHypo = new JLabel("Reconstruct hypothetical sequences:");
		lblRecHypo.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRecHypo.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRecHypo, "cell 4 20,alignx right");
		
		rdbtnRecHypoYes = new JRadioButton("Yes");
		rdbtnRecHypoYes.setBackground(new Color(204, 255, 255));
		rdbtnRecHypoYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRecHypoYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RecHypoToggle(false);
			}
		});
		panel.add(rdbtnRecHypoYes, "cell 5 20");
		
		rdbtnRecHypoNo = new JRadioButton("No");
		rdbtnRecHypoNo.setBackground(new Color(204, 255, 255));
		rdbtnRecHypoNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRecHypoNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RecHypoToggle(true);
			}
		});
		panel.add(rdbtnRecHypoNo, "cell 5 20");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInd, "cell 4 21,alignx right");
		
		rdbtnPrintIndYes = new JRadioButton("Yes");
		rdbtnPrintIndYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		panel.add(rdbtnPrintIndYes, "cell 5 21");
		
		rdbtnPrintIndNo = new JRadioButton("No");
		rdbtnPrintIndNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		panel.add(rdbtnPrintIndNo, "cell 5 21");
		
		btnDefaults = new JButton("Restore Defaults");
		btnDefaults.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				resetDefaults();
			}
		});
		btnDefaults.setFont(new Font("Arial", Font.BOLD, 13));	
		panel.add(btnDefaults, "cell 0 23");
		
		btnStored = new JButton("Read Init file");
		btnStored.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				getStoredSettings();
			}
		});
		btnStored.setFont(new Font("Arial", Font.BOLD, 13));	
		panel.add(btnStored, "cell 1 23");
		
		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				inputvals = getInputVals();
								
				btnExecute.setEnabled(false);	
				String title = "Dnaml Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread dnaMLThread = new Thread() {
						public void run() {
							runDnaMLThreads();
						}
			  	    };
			  	    dnaMLThread.start();
				}
				btnExecute.setEnabled(true);
				if(bootstrapCall)
				{
					saveSettings();
					final String[] bootstrap = new String[] {"Phylip", "bootstrap"};
					try {
						ConsenseUserInterface.main(bootstrap);
					} catch (Exception err) {
						err.printStackTrace();
					}						
					frmDnaMLControls.dispose();
				}
				else
				{
					btnDisplayTree.setEnabled(true);
					cmbxDisplayKind.setEnabled(true);
				}
			}
		});
		
		btnDisplayTree = new JButton("Display Tree In:");
		btnDisplayTree.setEnabled(false);
		btnDisplayTree.setHorizontalAlignment(SwingConstants.RIGHT);
		btnDisplayTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnDisplayTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				CallPlotType(cmbxDisplayKind.getSelectedIndex());
			}
		});
		panel.add(btnDisplayTree, "cell 2 23");

		cmbxDisplayKind = new JComboBox();
		cmbxDisplayKind.setEnabled(false);
		cmbxDisplayKind.setModel(new DefaultComboBoxModel(new String[] {"Drawgram", "Drawtree"}));
		cmbxDisplayKind.setSelectedIndex(0);
		cmbxDisplayKind.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxDisplayKind, "cell 3 23,growx");
		
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 5 23");
		
		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				saveSettings();
				if(phylipCall)
				{
					frmDnaMLControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 5 23");
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
		
		if (inputvals.InvarFract < 0.0) {
			String msg1 = "Input value: Fraction of sites invariant cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.InvarFract > 1.0) {
			String msg1 = "Input value: Fraction of sites invariant cannot greater than 1.0.";
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
	
	protected void runDnaMLThreads() {
    	try
    	{
       		// see if library exists
    		Native.loadLibrary("dnaml", DnaML.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("DnaML");
    		return;
    	}
		try 
		{
	  	    Thread dnaMLRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			DnaML dnaml = (DnaML) Native.loadLibrary("dnaml", DnaML.class);
		  	        dnaml.dnaml(		
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
		  	        		inputvals.SpeedAn,
		  	        		inputvals.GlobalRe,
		  	        		inputvals.RandInput,
		  	        		inputvals.RandNum,
		  	        		inputvals.Njumble,
		  	        		inputvals.OutRoot,
		  	        		inputvals.OutNum,
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
	  	    dnaMLRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (dnaMLRunThread.isAlive());
	  	    }
		} 
		catch (InterruptedException e) 
		{
			if (inputvals.PrintInd)
			{
				updateProgress();
			}
		}
	}
	
	protected DnaMLData getInputVals()
	{
		DnaMLData inputvals = new DnaMLData();
		// catch last array data entries
		SiteRateValues[Integer.parseInt(cmbxRateCatnum.getSelectedItem().toString())-1] = Double.parseDouble(txtSiteRate.getText());
		HMMRateValues[Integer.parseInt(cmbxHMMCat.getSelectedItem().toString())-1]  = Double.parseDouble(txtHMMRate.getText());
		HMMProbValues[Integer.parseInt(cmbxHMMCat.getSelectedItem().toString())-1]  = Double.parseDouble(txtHMMProb.getText());

		inputvals.infile = txtInputFile.getText();
		inputvals.intree = txtInputTree.getText();
		inputvals.wgtsfile = txtWeightFile.getText();
		inputvals.catsfile = txtCatFile.getText();
		inputvals.outfile = txtOutputFile.getText();
		inputvals.outfileopt = "w";
		inputvals.outtree = txtOutTree.getText();
		inputvals.outtreeopt = "w";
		if (cmbxTreeSearchMethod.getSelectedIndex() == 0)
		{
			inputvals.TreeUseMethod = "Yes";
		}
		else if (cmbxTreeSearchMethod.getSelectedIndex() == 1)
		{
			inputvals.TreeUseMethod = "No";
		}
		else //if (TreeSearchMethod.getSelectedIndex() == 2)
		{
			inputvals.TreeUseMethod = "rearrange";
		}
		inputvals.UseLengths = rdbtnUseLengthsYes.isSelected();
		inputvals.TTratio = Double.parseDouble(txtTTratio.getText());
		inputvals.useEmpBF = rdbtnEmpBFYes.isSelected();
		inputvals.BaseFreqA = Double.parseDouble(txtBaseFreqA.getText());
		inputvals.BaseFreqC = Double.parseDouble(txtBaseFreqC.getText());
		inputvals.BaseFreqG = Double.parseDouble(txtBaseFreqG.getText());
		inputvals.BaseFreqTU = Double.parseDouble(txtBaseFreqTU.getText());
		inputvals.OneCat = rdbtnCatYes.isSelected();
		inputvals.NumCats = Integer.parseInt(cmbxNumCat.getSelectedItem().toString());
		
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
		inputvals.SpeedAn = rdbtnSpeedAnNo.isSelected();
		inputvals.GlobalRe = rdbtnGlobalReYes.isSelected();
		inputvals.RandInput = rdbtnRandOrderYes.isSelected();
		inputvals.RandNum = Integer.parseInt(txtRandSeed.getText());
		inputvals.Njumble = Integer.parseInt(txtNumberJumble.getText());
		inputvals.OutRoot = rdbtnOutRootYes.isSelected();
		inputvals.OutNum = Integer.parseInt(txtOutRootNum.getText());
		inputvals.MultData = rdbtnAnalyzeMultDataYes.isSelected();
		inputvals.MultDSet = rdbtnDataSets.isSelected();
		inputvals.NumSeqs = Integer.parseInt(txtNumSeqs.getText());
		inputvals.InputSeq = rdbtnInputFileSeqYes.isSelected();
		inputvals.PrintData = rdbtnPrintDataYes.isSelected();
		inputvals.PrintInd = rdbtnPrintIndYes.isSelected();
		inputvals.PrintTree = rdbtnPrintTreeYes.isSelected();
		inputvals.WriteTree = rdbtnWriteTreeYes.isSelected();
		inputvals.DotDiff = rdbtnDotDiffYes.isSelected();
		inputvals.RecHypo = rdbtnRecHypoYes.isSelected();
		
	return inputvals;
	}
	
	protected void saveSettings(){
		inputvals = getInputVals();
		// there must be a better way to format this output, but this works for the prototype JRM
        try {
            BufferedWriter output = new BufferedWriter(new FileWriter(initFile));
            output.write("infile : "+inputvals.infile+"\n");
      		output.write("intree : "+inputvals.intree+"\n");
      		output.write("wgtsfile : "+inputvals.wgtsfile+"\n");
      		output.write("catsfile : "+inputvals.catsfile+"\n");
      		output.write("outfile : "+inputvals.outfile+"\n");
      		//output.write("outfileopt : "+inputvals.outfileopt+"\n"); //makes no sense to save, can change between runs JRM
       		output.write("outtree : "+inputvals.outtree+"\n");
      		//output.write("outtreeopt : "+inputvals.outtreeopt+"\n"); //makes no sense to save, can change between runs JRM
      		output.write("TreeUseMethod : "+inputvals.TreeUseMethod+"\n");
      		output.write("UseLengths : "+String.format("%b",inputvals.UseLengths)+"\n");
      		output.write("TTratio : "+String.format("%f",inputvals.TTratio)+"\n");
      		output.write("useEmpBF : "+String.format("%b",inputvals.useEmpBF)+"\n");
      		output.write("BaseFreqA : "+String.format("%f",inputvals.BaseFreqA)+"\n");
      		output.write("BaseFreqC : "+String.format("%f",inputvals.BaseFreqC)+"\n");
      		output.write("BaseFreqG : "+String.format("%f",inputvals.BaseFreqG)+"\n");
      		output.write("BaseFreqTU : "+String.format("%f",inputvals.BaseFreqTU)+"\n");
      		output.write("OneCat : "+String.format("%b",inputvals.OneCat)+"\n");
      		output.write("NumCats : "+inputvals.NumCats+"\n");
      		output.write("SiteRate1 : "+String.format("%f",inputvals.SiteRate1)+"\n");
            output.write("SiteRate2 : "+String.format("%f",inputvals.SiteRate2)+"\n");
            output.write("SiteRate3 : "+String.format("%f",inputvals.SiteRate3)+"\n");
            output.write("SiteRate4 : "+String.format("%f",inputvals.SiteRate4)+"\n");
            output.write("SiteRate5 : "+String.format("%f",inputvals.SiteRate5)+"\n");
            output.write("SiteRate6 : "+String.format("%f",inputvals.SiteRate6)+"\n");
            output.write("SiteRate7 : "+String.format("%f",inputvals.SiteRate7)+"\n");
            output.write("SiteRate8 : "+String.format("%f",inputvals.SiteRate8)+"\n");
            output.write("SiteRate9 : "+String.format("%f",inputvals.SiteRate9)+"\n");
      		output.write("RateVar : "+inputvals.RateVar+"\n");
      		output.write("AdjCor : "+String.format("%b",inputvals.AdjCor)+"\n");
      		output.write("BlockLen : "+String.format("%f",inputvals.BlockLen)+"\n");
      		output.write("CoeffVar : "+String.format("%f",inputvals.CoeffVar)+"\n");
      		output.write("NumRates : "+inputvals.NumRates+"\n");
      		output.write("HMMrate1 : "+String.format("%f",inputvals.HMMrate1)+"\n");
            output.write("HMMrate2 : "+String.format("%f",inputvals.HMMrate2)+"\n");
            output.write("HMMrate3 : "+String.format("%f",inputvals.HMMrate3)+"\n");
            output.write("HMMrate4 : "+String.format("%f",inputvals.HMMrate4)+"\n");
            output.write("HMMrate5 : "+String.format("%f",inputvals.HMMrate5)+"\n");
            output.write("HMMrate6 : "+String.format("%f",inputvals.HMMrate6)+"\n");
            output.write("HMMrate7 : "+String.format("%f",inputvals.HMMrate7)+"\n");
            output.write("HMMrate8 : "+String.format("%f",inputvals.HMMrate8)+"\n");
            output.write("HMMrate9 : "+String.format("%f",inputvals.HMMrate9)+"\n");
      		output.write("HMMprob1 : "+String.format("%f",inputvals.HMMprob1)+"\n");
            output.write("HMMprob2 : "+String.format("%f",inputvals.HMMprob2)+"\n");
            output.write("HMMprob3 : "+String.format("%f",inputvals.HMMprob3)+"\n");
            output.write("HMMprob4 : "+String.format("%f",inputvals.HMMprob4)+"\n");
            output.write("HMMprob5 : "+String.format("%f",inputvals.HMMprob5)+"\n");
            output.write("HMMprob6 : "+String.format("%f",inputvals.HMMprob6)+"\n");
            output.write("HMMprob7 : "+String.format("%f",inputvals.HMMprob7)+"\n");
            output.write("HMMprob8 : "+String.format("%f",inputvals.HMMprob8)+"\n");
            output.write("HMMprob9 : "+String.format("%f",inputvals.HMMprob9)+"\n");
            output.write("InvarFract : "+String.format("%f",inputvals.InvarFract)+"\n");
          	output.write("SitesWeight : "+String.format("%b",inputvals.SitesWeight)+"\n");
      		output.write("SpeedAn : "+String.format("%b",inputvals.SpeedAn)+"\n");
      		output.write("GlobalRe : "+String.format("%b",inputvals.GlobalRe)+"\n");
      		output.write("RandInput : "+String.format("%b",inputvals.RandInput)+"\n");
      		output.write("RandNum : "+inputvals.RandNum+"\n");
      		output.write("Njumble : "+inputvals.Njumble+"\n");
      		output.write("OutRoot : "+String.format("%b",inputvals.OutRoot)+"\n");
      		output.write("OutNum : "+inputvals.OutNum+"\n");
      		output.write("MultData : "+String.format("%b",inputvals.MultData)+"\n");
      		output.write("MultDSet : "+String.format("%b",inputvals.MultDSet)+"\n");
      		output.write("NumSeqs : "+inputvals.NumSeqs+"\n");
      		output.write("InputSeq : "+String.format("%b",inputvals.InputSeq)+"\n");
      		output.write("PrintData : "+String.format("%b",inputvals.PrintData)+"\n");
      		output.write("PrintInd : "+String.format("%b",inputvals.PrintInd)+"\n");
      		output.write("PrintTree : "+String.format("%b",inputvals.PrintTree)+"\n");
      		output.write("WriteTree : "+String.format("%b",inputvals.WriteTree)+"\n");
      		output.write("DotDiff : "+String.format("%b",inputvals.DotDiff)+"\n");
      		output.write("RecHypo : "+String.format("%b",inputvals.RecHypo)+"\n"); 
      		output.write("DisplayKind : "+cmbxDisplayKind.getSelectedIndex()+"\n");

            output.close();
        } catch ( IOException ioerr ) {
             ioerr.printStackTrace();
        }        
	}
	
	protected void getStoredSettings(){
		// because we are setting screen values directly, this is a tedious mess to set up JRM
	    try 
	    {
	    	Scanner scanner =  new Scanner(new File(initFile));
	        while (scanner.hasNextLine()){
	        	Scanner linescan =  new Scanner( scanner.nextLine());
	        	linescan.useDelimiter(" : ");
	        	String label = linescan.next();
	        	String value = linescan.next();
	        	if ("infile".equals(label)){
	        		txtInputFile.setText(value);
	        	}
	     		else if ("intree".equals(label)){
	     			txtInputTree.setText(value);
	     		}
	    		else if ("wgtsfile".equals(label)){
	    			txtWeightFile.setText(value);
	    		}
				else if ("catsfile".equals(label)){
					txtCatFile.setText(value);
				}
				else if ("outfile".equals(label)){
					txtOutputFile.setText(value);
				}
				else if ("outtree".equals(label)){
					txtOutTree.setText(value);
				}
				else if ("TreeUseMethod".equals(label)){
					if (value.contains("rearrange"))
					{
						cmbxTreeSearchMethod.setSelectedIndex(2);
					}
					else if (value.contains("No"))
					{
						cmbxTreeSearchMethod.setSelectedIndex(1);					
					}
					else// if (value.contains("Yes"))
					{
						cmbxTreeSearchMethod.setSelectedIndex(0);					
					}
				}
	    		else if ("UseLengths".equals(label)){
	    			if ("true".equals(value))
	    			{
	     				rdbtnUseLengthsYes.setSelected(true);
	     				rdbtnUseLengthsNo.setSelected(false);
	     			}
	    			else
	    			{
	     				rdbtnUseLengthsYes.setSelected(false);
	     				rdbtnUseLengthsNo.setSelected(true);    				
	    			}
	    		}
	    		else if ("TTratio".equals(label)){
	    			txtTTratio.setText(value);
	    		}
	    		else if ("useEmpBF".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				EmpBFToggle(true);
	    			}
	    			else
	    			{    				
	    				EmpBFToggle(false);
	    			}
	    		}
	    		else if ("BaseFreqA".equals(label)){
	    			txtBaseFreqA.setText(value);
	    		}
	    		else if ("BaseFreqC".equals(label)){
	    			txtBaseFreqC.setText(value);
	    		}
	    		else if ("BaseFreqG".equals(label)){
	    			txtBaseFreqG.setText(value);
	    		}
	    		else if ("BaseFreqTU".equals(label)){
	    			txtBaseFreqTU.setText(value);
	    		}
	    		else if ("OneCat".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				CatToggle(true);
	    			}
	    			else
	    			{    				
	    				CatToggle(false);
	    			}
	    		}
	    		else if ("NumCats".equals(label)){
	    			cmbxNumCat.setSelectedItem(Integer.parseInt(value)-1);
	    		}
	    		else if ("SiteRate1".equals(label)){
	    			SiteRateValues[0] = Double.parseDouble(value);
	    		}    
	    		else if ("SiteRate2".equals(label)){
	    			SiteRateValues[1] = Double.parseDouble(value);
	    		}    
	    		else if ("SiteRate3".equals(label)){
	    			SiteRateValues[2] = Double.parseDouble(value);
	    		}    
	    		else if ("SiteRate4".equals(label)){
	    			SiteRateValues[3] = Double.parseDouble(value);
	    		}    
	    		else if ("SiteRate5".equals(label)){
	    			SiteRateValues[4] = Double.parseDouble(value);
	    		}    
	    		else if ("SiteRate6".equals(label)){
	    			SiteRateValues[5] = Double.parseDouble(value);
	    		}    
	    		else if ("SiteRate7".equals(label)){
	    			SiteRateValues[6] = Double.parseDouble(value);
	    		}    
	    		else if ("SiteRate8".equals(label)){
	    			SiteRateValues[8] = Double.parseDouble(value);
	    		}   
	    		else if ("SiteRate9".equals(label)){
	    			SiteRateValues[8] = Double.parseDouble(value);
	    		}  
	    		else if ("RateVar".equals(label)){
	    			if("constant".equals(value))
	    			{
	    				cmbxRateSite.setSelectedIndex(0);
	    			}
	    			else if("gamma".equals(value))
	    			{
	    				cmbxRateSite.setSelectedIndex(1);
	    			}
	    			else if("invar".equals(value))
	    			{
	    				cmbxRateSite.setSelectedIndex(2);
	    			}
	    			else //if("hmm".equals(value))
	    			{
	    				cmbxRateSite.setSelectedIndex(3);
	    			}
	    		}
	    		else if ("AdjCor".equals(label)){
	    			if ("true".equals(value))
	    			{
	       			 	rdbtnAdjCorYes.setSelected(true);
	       			 	rdbtnAdjCorNo.setSelected(false);
	     			}
	    			else
	    			{    				
	       			 	rdbtnAdjCorYes.setSelected(false);
	       			 	rdbtnAdjCorNo.setSelected(true);
	    			}
	    		}
	    		else if ("BlockLen".equals(label)){
	    			txtBlockLen.setText(value);
	    		}
	    		else if ("CoeffVar".equals(label)){
	    			txtCoeffVar.setText(value);
	    		}
	    		else if ("NumRates".equals(label)){
	    			cmbxHMMcount.setSelectedItem(Integer.parseInt(value)-1);
	    		}
	    		else if ("HMMrate1".equals(label)){
	            	HMMRateValues[0] = Double.parseDouble(value);
	            	}
	            else if ("HMMrate2".equals(label)){
	            	HMMRateValues[1] = Double.parseDouble(value);
	            	}
	            else if ("HMMrate3".equals(label)){
	            	HMMRateValues[2] = Double.parseDouble(value);
	            	}
	            else if ("HMMrate4".equals(label)){
	            	HMMRateValues[3] = Double.parseDouble(value);
	            	}
	            else if ("HMMrate5".equals(label)){
	            	HMMRateValues[4] = Double.parseDouble(value);
	            	}
	            else if ("HMMrate6".equals(label)){
	            	HMMRateValues[5] = Double.parseDouble(value);
	            	}
	            else if ("HMMrate7".equals(label)){
	            	HMMRateValues[6] = Double.parseDouble(value);
	            	}
	            else if ("HMMrate8".equals(label)){
	            	HMMRateValues[7] = Double.parseDouble(value);
	            	}
	            else if ("HMMrate9".equals(label)){
	            	HMMRateValues[8] = Double.parseDouble(value);
	            	}
	    		else if ("HMMprob1".equals(label)){
	            	HMMProbValues[0] = Double.parseDouble(value);
	            	}
	            else if ("HMMprob2".equals(label)){
	            	HMMProbValues[1] = Double.parseDouble(value);
	            	}
	            else if ("HMMprob3".equals(label)){
	            	HMMProbValues[2] = Double.parseDouble(value);
	            	}
	            else if ("HMMprob4".equals(label)){
	            	HMMProbValues[3] = Double.parseDouble(value);
	            	}
	            else if ("HMMprob5".equals(label)){
	            	HMMProbValues[4] = Double.parseDouble(value);
	            	}
	            else if ("HMMprob6".equals(label)){
	            	HMMProbValues[5] = Double.parseDouble(value);
	            	}
	            else if ("HMMprob7".equals(label)){
	            	HMMProbValues[6] = Double.parseDouble(value);
	            	}
	            else if ("HMMprob8".equals(label)){
	             	HMMProbValues[7] = Double.parseDouble(value);
	            }
	            else if ("HMMprob9".equals(label)){
	            	HMMProbValues[8] = Double.parseDouble(value);
	            }
	            else if ("InvarFract".equals(label)){
	            	txtFracInvar.setText(value);
	            }
	    		else if ("SitesWeight".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				SitesWeightToggle(false);
	    			}
	    			else
	    			{    				
	    				SitesWeightToggle(true);
	    			}
	    		}
	    		else if ("SpeedAn".equals(label)){
	    			if ("true".equals(value))
	    			{
	       				rdbtnSpeedAnNo.setSelected(true);
	       				rdbtnSpeedAnYes.setSelected(false);
	    			}
	    			else
	    			{    				
	       				rdbtnSpeedAnNo.setSelected(false);
	       				rdbtnSpeedAnYes.setSelected(true);
	    			}
	    		}
	    		else if ("GlobalRe".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnGlobalReYes.setSelected(true);
	    				rdbtnGlobalReNo.setSelected(false);
	    			}
	    			else
	    			{    				
	    				rdbtnGlobalReYes.setSelected(false);
	    				rdbtnGlobalReNo.setSelected(true);
	    			}
	    		}
	    		else if ("RandInput".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnRandOrderYes.setSelected(true);
	    				rdbtnRandOrderNo.setSelected(false);
	    			}
	    			else
	    			{    				
	    				rdbtnRandOrderYes.setSelected(false);
	    				rdbtnRandOrderNo.setSelected(true);
	    			}
	    		}
	    		else if ("RandNum".equals(label)){
	    			txtRandSeed.setText(value);
	    		}
	    		else if ("Njumble".equals(label)){
	    			txtNumberJumble.setText(value);
	    		}
	    		else if ("OutRoot".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				OutRootToggle(false);
	    			}
	    			else
	    			{    				
	    				OutRootToggle(true);
	    			}
	    		}
	    		else if ("OutNum".equals(label)){
	    			txtOutRootNum.setText(value);
	    		}
	    		else if ("MultData".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				MultToggle(true);	
	    			}
	    			else
	    			{    				
	    				MultToggle(false);	    			
	    			}
	    		}
	    		else if ("MultDSet".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnDataSets.setSelected(true);
	     			}
	    			else
	    			{    				
	       				rdbtnDataSets.setSelected(false);
	    			}
	    		}
	    		else if ("NumSeqs".equals(label)){
	    			txtNumSeqs.setText(value);
	    		}
	    		else if ("InputSeq".equals(label)){
	    			if ("true".equals(value))
	    			{
	       				rdbtnInputFileSeqYes.setSelected(true);
	       				rdbtnInputFileSeqNo.setSelected(false);
	       			}
	    			else
	    			{    				
	      				rdbtnInputFileSeqYes.setSelected(false);
	       				rdbtnInputFileSeqNo.setSelected(true);
	    			}
	    		}
	    		else if ("PrintData".equals(label)){
	    			if ("true".equals(value))
	    			{
	       				rdbtnPrintDataYes.setSelected(true);
	       				rdbtnPrintDataNo.setSelected(false);
	       			}
	    			else
	    			{    				
	       				rdbtnPrintDataYes.setSelected(false);
	       				rdbtnPrintDataNo.setSelected(true);
	    			}
	    		}
	    		else if ("PrintInd".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnPrintIndYes.setSelected(true);
	    				rdbtnPrintIndNo.setSelected(false);
	    			}
	    			else
	    			{    				
	    				rdbtnPrintIndYes.setSelected(false);
	    				rdbtnPrintIndNo.setSelected(true);
	    			}
	    		}
	    		else if ("PrintTree".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnPrintTreeYes.setSelected(true);
	    				rdbtnPrintTreeNo.setSelected(false);
	    			}
	    			else
	    			{    				
	    				rdbtnPrintTreeYes.setSelected(false);
	    				rdbtnPrintTreeNo.setSelected(true);
	    			}
	    		}
	    		else if ("DotDiff".equals(label)){
	    			if ("true".equals(value))
	    			{
	      				rdbtnDotDiffYes.setSelected(true);
	       				rdbtnDotDiffNo.setSelected(false);
	     			}
	    			else
	    			{    				
	       				rdbtnDotDiffYes.setSelected(false);
	       				rdbtnDotDiffNo.setSelected(true);
	    			}
	    		}
	    		else if ("WriteTree".equals(label)){
	    			if ("true".equals(value))
	    			{
	       				rdbtnWriteTreeYes.setSelected(true);
	    				rdbtnWriteTreeNo.setSelected(false);
	      			}
	    			else
	    			{    				
	    				rdbtnWriteTreeYes.setSelected(false);
	    				rdbtnWriteTreeNo.setSelected(true);
	    			}
	
	    		}
	    		else if ("RecHypo".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnRecHypoYes.setSelected(true);
	    				rdbtnRecHypoNo.setSelected(false);
	    			}
	    			else
	    			{    				
	       				rdbtnRecHypoYes.setSelected(false);
	    				rdbtnRecHypoNo.setSelected(true);
	    			}
	    		}
	    		else if ("DisplayKind".equals(label)){
	    			cmbxDisplayKind.setSelectedIndex(Integer.parseInt(value));
	    		}
	    		else {
	    			String msg = "Unknown label: ";
	    			msg += label;
	    			msg += " with value: ";
	    			msg += value;
	    			msg += " found in dnamlInit.txt.";
	    			JOptionPane.showMessageDialog(null, msg, "Warning", JOptionPane.WARNING_MESSAGE);			
	    		}
	        }  
	        // set to default values
	        txtSiteRate.setText(Double.toString(SiteRateValues[0]));
			txtHMMRate.setText(Double.toString(HMMRateValues[0]));
			txtHMMProb.setText(Double.toString(HMMProbValues[0]));

	    }
		catch (FileNotFoundException e)
		{
			// if it's not there, use the defaults
			resetDefaults();
		} 	
	}
	
	protected void resetDefaults()
	{
		// reset DnaML to default values
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
		
		HMMProbValues[0] = 1.0;
		HMMProbValues[1] = 1.0;
		HMMProbValues[2] = 1.0;
		HMMProbValues[3] = 1.0;
		HMMProbValues[4] = 1.0;
		HMMProbValues[5] = 1.0;
		HMMProbValues[6] = 1.0;
		HMMProbValues[7] = 1.0;
		HMMProbValues[8] = 1.0;
		
		txtInputFile.setText("infile");
		btnInputTree.setEnabled(false);
		txtInputTree.setEnabled(false);
		txtInputTree.setText("intree");
		btnWeightFile.setEnabled(false);
		txtWeightFile.setEnabled(false);
		txtWeightFile.setText("weightfile");
		btnCatFile.setEnabled(false);
		txtCatFile.setEnabled(false);
		txtCatFile.setText("catfile");
		txtOutputFile.setText("outfile");
		txtOutTree.setText("outtree");
		cmbxTreeSearchMethod.setSelectedIndex(0);
		lblUseLengths.setEnabled(false);
		rdbtnUseLengthsYes.setEnabled(false);
		rdbtnUseLengthsYes.setSelected(false);
		rdbtnUseLengthsNo.setSelected(true);
		rdbtnUseLengthsNo.setEnabled(false);
		txtTTratio.setText("2.00");
		rdbtnEmpBFYes.setSelected(true);
		rdbtnEmpBFNo.setSelected(false);
		lblBaseFreq.setEnabled(false);
		lblBaseFreqA.setEnabled(false);
		txtBaseFreqA.setEnabled(false);
		txtBaseFreqA.setText("0.25");
		lblBaseFreqC.setEnabled(false);
		txtBaseFreqC.setEnabled(false);
		txtBaseFreqC.setText("0.25");
		lblBaseFreqG.setEnabled(false);
		txtBaseFreqG.setEnabled(false);
		txtBaseFreqG.setText("0.25");
		lblBaseFreqTU.setEnabled(false);
		txtBaseFreqTU.setEnabled(false);
		txtBaseFreqTU.setText("0.25");
		rdbtnCatYes.setSelected(true);
		rdbtnCatNo.setSelected(false);
		lblNumCat.setEnabled(false);
		cmbxNumCat.setEnabled(false);
		cmbxNumCat.setSelectedIndex(0);
		lblcatnum.setEnabled(false);
		cmbxRateCatnum.setEnabled(false);
		lblRate.setEnabled(false);
		txtSiteRate.setText("1.0");
		txtSiteRate.setEnabled(false);
		cmbxRateSite.setSelectedIndex(0);
		lblRatesAdjCor.setEnabled(false);
		rdbtnAdjCorYes.setEnabled(false);
		rdbtnAdjCorYes.setSelected(false);
		rdbtnAdjCorNo.setSelected(true);
		rdbtnAdjCorNo.setEnabled(false);
		lblMeanLen.setEnabled(false);
		txtBlockLen.setText("1");
		txtBlockLen.setEnabled(false);
		lblCoeffVar.setEnabled(false);
		txtCoeffVar.setEnabled(false);
		txtCoeffVar.setText("1.0");
		lblCoeffVarNote.setVisible(false);
		lblHMM.setEnabled(false);
		lblInvar.setVisible(false);
		lblNumHMMcats.setEnabled(false);
		cmbxHMMcount.setSelectedIndex(0);
		cmbxHMMcount.setEnabled(false);
		lblHMMCat.setEnabled(false);
		cmbxHMMCat.setEnabled(false);
		lblHMMRate.setEnabled(false);
		txtHMMRate.setText("1.0");
		txtHMMRate.setEnabled(false);
		lblProb.setEnabled(false);
		txtHMMProb.setText("1.0");
		txtHMMProb.setEnabled(false);
		lblFracInvar.setEnabled(false);
		txtFracInvar.setText("1.0");
		txtFracInvar.setEnabled(false);	
		rdbtnSitesWeightYes.setSelected(false);
		rdbtnSitesWeightNo.setSelected(true);
		rdbtnSpeedAnYes.setSelected(true);
		rdbtnSpeedAnNo.setSelected(false);
		rdbtnGlobalReYes.setSelected(false);
		rdbtnGlobalReNo.setSelected(true);
		rdbtnRandOrderYes.setSelected(false);
		rdbtnRandOrderNo.setSelected(true);
		lblRandSeed.setEnabled(false);
		txtRandSeed.setText("1");
		txtRandSeed.setEnabled(false);
		lblRandOdd.setEnabled(false);
		lblNumberJumble.setEnabled(false);
		txtNumberJumble.setText("1");
		txtNumberJumble.setEnabled(false);
		rdbtnOutRootYes.setSelected(false);
		rdbtnOutRootNo.setSelected(true);
		lblOutRootNum.setEnabled(false);
		txtOutRootNum.setText("1");
		txtOutRootNum.setEnabled(false);
		rdbtnAnalyzeMultDataYes.setSelected(false);
		rdbtnAnalyzeMultDataNo.setSelected(true);
		lblMultData.setEnabled(false);
		rdbtnDataSets.setSelected(true);
		rdbtnDataSets.setEnabled(false);
		rdbtnWeights.setSelected(false);
		rdbtnWeights.setEnabled(false);
		lblHowManyData.setEnabled(false);
		txtNumSeqs.setText("1");
		txtNumSeqs.setEnabled(false);
		lblInputSeq.setEnabled(false);
		rdbtnInputFileSeqYes.setEnabled(false);
		rdbtnInputFileSeqYes.setSelected(true);
		rdbtnInputFileSeqNo.setEnabled(false);
		rdbtnInputFileSeqNo.setSelected(false);
		rdbtnPrintDataYes.setSelected(false);
		rdbtnPrintDataNo.setSelected(true);
		rdbtnPrintTreeYes.setSelected(true);
		rdbtnPrintTreeNo.setSelected(false);
		rdbtnWriteTreeYes.setSelected(true);
		rdbtnWriteTreeNo.setSelected(false);
		rdbtnDotDiffYes.setSelected(true);
		rdbtnDotDiffNo.setSelected(false);
		rdbtnRecHypoYes.setSelected(false);
		rdbtnRecHypoNo.setSelected(true);
		rdbtnPrintIndYes.setSelected(true);
		rdbtnPrintIndNo.setSelected(false);
		btnDisplayTree.setEnabled(false);
		cmbxDisplayKind.setEnabled(false);

	}
	  
	protected void getSeqBootSettings()
	{
		// get the parameters needed
	    try 
	    {
	    	Scanner scanner =  new Scanner(new File("seqbootInit.txt"));
	        while (scanner.hasNextLine()){
	        	Scanner linescan =  new Scanner( scanner.nextLine());
	        	linescan.useDelimiter(" : ");
	        	String label = linescan.next();
	        	String value = linescan.next();
	      		if("outfile".equals(label)){
	      			txtInputFile.setText(value);
	        	}
	      		else if("Replicates".equals(label)){
	      			txtNumSeqs.setText(value);
	        	}
	        }
	    }
		catch (FileNotFoundException e)
		{
			String msg = "Input file: seqbootInit.txt does not exist.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);			
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