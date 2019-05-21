package phylip;
import java.awt.EventQueue;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import com.sun.jna.Library;
import com.sun.jna.Native;

import java.awt.Font;
import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class CodMLUserInterface {
   public interface CodML extends Library {
        public void codml(
        		String infile,
        		String intree,
        		String wgtsfile,
        		String outfile,
        		String outfileopt,
        		String outtree,
        		String outtreeopt,
        		String TreeUseMethod,
        		boolean UseLengths,
        		boolean InputNuc,
        		String NucSubModel,
        		double TTratio,
        		double NSratio,
        		boolean useEmpBF,
        		double BaseFreqA,
        		double BaseFreqC,
        		double BaseFreqG,
        		double BaseFreqTU,
        		String GeneCode,
        		boolean OneOmega,
        		boolean AdjOmegasCor,
        		int OmegaBlockLen,
        		boolean CodonsWgted,
        		int NumOmegas,
        		
        		// these are explicitly named because JNA doesn't pass arrays gracefully
        		double OmegaVal1,
                double OmegaVal2,
                double OmegaVal3,
                double OmegaVal4,
                double OmegaVal5,
                double OmegaVal6,
                double OmegaVal7,
                double OmegaVal8,
                double OmegaVal9,
                //
        		// these are explicitly named because JNA doesn't pass arrays gracefully
                double OmegaProb1,
                double OmegaProb2,
                double OmegaProb3,
                double OmegaProb4,
                double OmegaProb5,
                double OmegaProb6,
                double OmegaProb7,
                double OmegaProb8,
                double OmegaProb9,
                
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
   public class CodMLData {
		String infile;
		String intree;
		String wgtsfile;
		String outfile;
		String outfileopt;
		String outtree;
		String outtreeopt;
		String TreeUseMethod;
		boolean UseLengths;
		boolean InputNuc;
		String NucSubModel;
		double TTratio;
		double NSratio;
		boolean useEmpBF;
		double BaseFreqA;
		double BaseFreqC;
		double BaseFreqG;
		double BaseFreqTU;
		String GeneCode;
		boolean OneOmega;
		boolean AdjOmegasCor;
		int OmegaBlockLen;
		boolean CodonsWgted;
		int NumOmegas;
		
		// these are explicitly named because JNA doesn't pass arrays gracefully
		double OmegaVal1;
		double OmegaVal2;
		double OmegaVal3;
		double OmegaVal4;
		double OmegaVal5;
		double OmegaVal6;
		double OmegaVal7;
		double OmegaVal8;
		double OmegaVal9;

		// these are explicitly named because JNA doesn't pass arrays gracefully
	    double OmegaProb1;
	    double OmegaProb2;
	    double OmegaProb3;
	    double OmegaProb4;
	    double OmegaProb5;
	    double OmegaProb6;
	    double OmegaProb7;
	    double OmegaProb8;
	    double OmegaProb9;
   
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
   
	private CodMLData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean ExplicitWgts;
	private double[] OmegaVal;
	private double[] OmegaProb;
	private int lasthmmcat;
	private boolean phylipCall;

	private JFrame frmCodMLControls;
	private JTextField txtInputFile;
	private JButton btnInputFile;
	private JTextField txtInputTree;
	private JButton btnInputTree;
	private JButton btnExecute;
	private JButton btnQuit;
	private JTextField txtOutputFile;
	private JButton btnOutputFile;
	private JButton btnWeightFile;
	private JTextField txtWeightFile;
	private JButton btnOutputTree;
	private JTextField txtOutputTree;
	private JLabel lblSearchBest;
	private JComboBox cmbxTreeSearchMethod;
	private JLabel lblUseLengths;
	private JRadioButton rdbtnUseLengthsYes;
	private JRadioButton rdbtnUseLengthsNo;
	private JLabel lblDotDiff;
	private JRadioButton rdbtnDotDiffYes;
	private JRadioButton rdbtnDotDiffNo;
	private JLabel lblTTratio;
	private JTextField txtTTratio;
	private JLabel lblNSratio;
	private JTextField txtNSratio;
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
	
	private JScrollPane scrollPane;
	private JPanel panel;
	private JLabel lblReadInputData;
	private JRadioButton rdbtnNucleotides;
	private JRadioButton rdbtnAminoAcids;
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
	private JLabel lblNucSubModel;
	private JComboBox cmbxNucSubMethod;
	private JLabel lblCodonTable;
	private JComboBox cmbxCodonTable;
	private JLabel lblOneOmega;
	private JRadioButton rdbtnOneOmegaYes;
	private JRadioButton rdbtnOneOmegaNo;
	private JLabel lblOmegasCor;
	private JRadioButton rdbtnOmegasCorYes;
	private JRadioButton rdbtnOmegasCorNo;
	private JLabel lblBlockLen;
	private JTextField txtBlockLen;
	private JLabel lblGt1;
	private JLabel lblCodonsWeighted;
	private JRadioButton rdbtnCodonsWeightedYes;
	private JRadioButton rdbtnCodonsWeightedNo;
	private JLabel lblSpeedyAnalysis;
	private JRadioButton rdbtnSpeedyYes;
	private JRadioButton rdbtnSpeedyNo;
	private JLabel lblGlobalRearr;
	private JRadioButton rdbtnGlobalRearrYes;
	private JRadioButton rdbtnGlobalRearrNo;
	private JLabel lblRandInSeqOrd;
	private JRadioButton rdbtnRandInSeqOrdYes;
	private JRadioButton rdbtnRandInSeqOrdNo;
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
	private JLabel lblRandOdd;
	private JLabel lblNumberJumble;
	private JTextField txtNumberJumble;
	private JLabel lblHowManyData;
	private JTextField txtNumSeqs;
	private JLabel lblOmegaValuesHMM;
	private JLabel lblNumCats;
	private JComboBox cmbxNumCats;
	private JLabel lblCategory;
	private JComboBox cmbxCategory;
	private JLabel lblRate;
	private JTextField txtOmegaVal;
	private JLabel lblProbability;
	private JTextField txtOmegaProb;

	/**
	 * Launch ttnNewButton;

	/**
	 * Launch the application
	 * 	public static void main(String[] args) {.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					CodMLUserInterface window = new CodMLUserInterface(args);
					window.frmCodMLControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmCodMLControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}
	
	protected void IntreeToggle(int selected){
		if(selected == 0) //yes
		{	
			txtInputTree.setEnabled(false);
			btnInputTree.setEnabled(false);
			lblUseLengths.setEnabled(false);
			rdbtnUseLengthsYes.setEnabled(false);
			rdbtnUseLengthsNo.setEnabled(false);
			lblRandInSeqOrd.setEnabled(true);
			rdbtnRandInSeqOrdYes.setEnabled(true);
			rdbtnRandInSeqOrdNo.setEnabled(true);	
			RandOrderToggle(rdbtnRandInSeqOrdYes.isSelected());
		}
		else // no or rearrange
		{	
			txtInputTree.setEnabled(true);
			btnInputTree.setEnabled(true);
			lblRandInSeqOrd.setEnabled(false);
			rdbtnRandInSeqOrdYes.setEnabled(false);
			rdbtnRandInSeqOrdNo.setEnabled(false);	
			lblRandSeed.setEnabled(false);
			txtRandSeed.setEnabled(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEnabled(false);
			lblRandOdd.setEnabled(false);
			if (selected == 1) // no
			{
				lblUseLengths.setEnabled(true);
				rdbtnUseLengthsYes.setEnabled(true);
				rdbtnUseLengthsNo.setEnabled(true);
			}
			else // rearrange
			{
				lblUseLengths.setEnabled(false);
				rdbtnUseLengthsYes.setEnabled(false);
				rdbtnUseLengthsNo.setEnabled(false);
			}
		}
	}
	
	protected void NucSubMethodToggle(int selected){
		if(selected == 0) //f84
		{				
			lblTTratio.setEnabled(true);
			txtTTratio.setEnabled(true);
			lblEmpBF.setEnabled(true);
			rdbtnEmpBFYes.setEnabled(true);
			rdbtnEmpBFNo.setEnabled(true);
			if (rdbtnEmpBFYes.isSelected())
			{
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
			else
			{
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
		else if(selected == 1) //HKY
		{	
			lblTTratio.setEnabled(true);
			txtTTratio.setEnabled(true);
			lblEmpBF.setEnabled(true);
			rdbtnEmpBFYes.setEnabled(true);
			rdbtnEmpBFNo.setEnabled(true);
			if (rdbtnEmpBFYes.isSelected())
			{
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
			else
			{
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
		else if(selected == 2) //Kimura 2-parameter
		{	
			lblTTratio.setEnabled(true);
			txtTTratio.setEnabled(true);
			lblEmpBF.setEnabled(false);
			rdbtnEmpBFYes.setEnabled(false);
			rdbtnEmpBFNo.setEnabled(false);
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
		else // if(selected == 3) //Jukes-Cantor
		{	
			lblTTratio.setEnabled(false);
			txtTTratio.setEnabled(false);
			lblEmpBF.setEnabled(false);
			rdbtnEmpBFYes.setEnabled(false);
			rdbtnEmpBFNo.setEnabled(false);
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
			//txtOutTree.setEnabled(true);
		}
		else{
			rdbtnWriteTreeYes.setSelected(false);
			rdbtnWriteTreeNo.setSelected(true);
			btnOutputTree.setEnabled(false);
			//txtOutTree.setEnabled(false);
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

	protected void DotDiffToggle(boolean isUseDot) {
		if (isUseDot) {
			rdbtnDotDiffYes.setSelected(true);
			rdbtnDotDiffNo.setSelected(false);
		} else {
			rdbtnDotDiffYes.setSelected(false);
			rdbtnDotDiffNo.setSelected(true);
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
			}
			else
			{
				btnWeightFile.setEnabled(false);
				txtWeightFile.setEnabled(false);	
			}
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
			}
		}
		else{
			rdbtnDataSets.setSelected(false);
			rdbtnWeights.setSelected(true);
			RandOrderToggle(false);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
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
	
	protected void RandOrderToggle(boolean isRand){
		if(isRand){
			rdbtnRandInSeqOrdNo.setSelected(false);
			rdbtnRandInSeqOrdYes.setSelected(true);
			lblRandSeed.setEnabled(true);
			txtRandSeed.setEnabled(true);
			lblNumberJumble.setEnabled(true);
			txtNumberJumble.setEnabled(true);
			lblRandOdd.setEnabled(true);
		}
		else{
			rdbtnRandInSeqOrdNo.setSelected(true);
			rdbtnRandInSeqOrdYes.setSelected(false);
			lblRandSeed.setEnabled(false);
			txtRandSeed.setEnabled(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEnabled(false);
			lblRandOdd.setEnabled(false);
		}
	}
	
	protected void NucleotidesToggle(boolean isNuc){
		if(isNuc){		
			rdbtnNucleotides.setSelected(true);
			rdbtnAminoAcids.setSelected(false);
		}
		else{
			rdbtnNucleotides.setSelected(false);
			rdbtnAminoAcids.setSelected(true);
		}
	}
	
	protected void OneOmegaToggle(boolean isYes){
		if(isYes){		
			rdbtnOneOmegaYes.setSelected(true);
			rdbtnOneOmegaNo.setSelected(false);
			lblOmegasCor.setEnabled(false);
			rdbtnOmegasCorYes.setEnabled(false);
			rdbtnOmegasCorNo.setEnabled(false);
			lblBlockLen.setEnabled(false);
			txtBlockLen.setEnabled(false);
			lblGt1.setEnabled(false);
			lblNSratio.setEnabled(true);
			txtNSratio.setEnabled(true);
		}
		else{
			rdbtnOneOmegaYes.setSelected(false);
			rdbtnOneOmegaNo.setSelected(true);
			lblOmegasCor.setEnabled(true);
			rdbtnOmegasCorYes.setEnabled(true);
			rdbtnOmegasCorNo.setEnabled(true);
			lblNSratio.setEnabled(false);
			txtNSratio.setEnabled(false);
			OmegaCorToggle(rdbtnOmegasCorYes.isSelected());				
		}
	}
	
	protected void OmegaCorToggle(boolean isCor)
	{
		if(isCor)
		{
			rdbtnOmegasCorYes.setSelected(true);
			rdbtnOmegasCorNo.setSelected(false);
			lblBlockLen.setEnabled(true);
			txtBlockLen.setEnabled(true);
			lblGt1.setEnabled(true);
		}
		else
		{
			rdbtnOmegasCorYes.setSelected(false);
			rdbtnOmegasCorNo.setSelected(true);
			lblBlockLen.setEnabled(false);
			txtBlockLen.setEnabled(false);
			lblGt1.setEnabled(false);
		}
	}
	
	protected void SpeedyToggle(boolean isSpeedy){
		if(isSpeedy){
			rdbtnSpeedyYes.setSelected(true);
			rdbtnSpeedyNo.setSelected(false);
		}
		else{
			rdbtnSpeedyYes.setSelected(false);
			rdbtnSpeedyNo.setSelected(true);
		}
	}
	
	protected void CodonWgtToggle(boolean isWgt){
		if(isWgt){
			rdbtnCodonsWeightedYes.setSelected(true);
			rdbtnCodonsWeightedNo.setSelected(false);
			lblOmegaValuesHMM.setEnabled(true);
			lblNumCats.setEnabled(true);
			cmbxNumCats.setEnabled(true);
			lblCategory.setEnabled(true);
			cmbxCategory.setEnabled(true);
			lblRate.setEnabled(true);
			txtOmegaVal.setEnabled(true);
			lblProbability.setEnabled(true);
			txtOmegaProb.setEnabled(true);
		}
		else{
			rdbtnCodonsWeightedYes.setSelected(false);
			rdbtnCodonsWeightedNo.setSelected(true);
			lblOmegaValuesHMM.setEnabled(false);
			lblNumCats.setEnabled(false);
			cmbxNumCats.setEnabled(false);
			lblCategory.setEnabled(false);
			cmbxCategory.setEnabled(false);
			lblRate.setEnabled(false);
			txtOmegaVal.setEnabled(false);
			lblProbability.setEnabled(false);
			txtOmegaProb.setEnabled(false);
		}
	}
	
	protected void GlobalRearrToggle(boolean isGlobal){
		if(isGlobal){
			rdbtnGlobalRearrYes.setSelected(true);
			rdbtnGlobalRearrNo.setSelected(false);
		}
		else{
			rdbtnGlobalRearrYes.setSelected(false);
			rdbtnGlobalRearrNo.setSelected(true);
		}
	}
	
	protected void DisplayHMMValues(int catnum){
		if ((catnum-1) != lasthmmcat){
			OmegaVal[lasthmmcat] = Double.parseDouble(txtOmegaVal.getText());
			String OmegaValstr = Double.toString(OmegaVal[catnum-1]);
			txtOmegaVal.setText(OmegaValstr);
			
			OmegaProb[lasthmmcat] = Double.parseDouble(txtOmegaProb.getText());
			String ratestr = Double.toString(OmegaProb[catnum-1]);
			txtOmegaProb.setText(ratestr);
		}
		lasthmmcat = catnum-1;
	}

	/**
	 * Create the application.
	 */
	public CodMLUserInterface(String[] args) {
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
		
		OmegaVal = new double[9];
		OmegaVal[0] = 1.0;
		OmegaVal[1] = 1.0;
		OmegaVal[2] = 1.0;
		OmegaVal[3] = 1.0;
		OmegaVal[4] = 1.0;
		OmegaVal[5] = 1.0;
		OmegaVal[6] = 1.0;
		OmegaVal[7] = 1.0;
		OmegaVal[8] = 1.0;
		
		lasthmmcat = 0;
		
		OmegaProb = new double[9];
		OmegaProb[0] = 1.0;
		OmegaProb[1] = 1.0;
		OmegaProb[2] = 1.0;
		OmegaProb[3] = 1.0;
		OmegaProb[4] = 1.0;
		OmegaProb[5] = 1.0;
		OmegaProb[6] = 1.0;
		OmegaProb[7] = 1.0;
		OmegaProb[8] = 1.0;
	
		filedir = System.getProperty("user.dir");


		frmCodMLControls = new JFrame();
		frmCodMLControls.setBackground(new Color(204, 255, 255));
		frmCodMLControls.setTitle("Codml");
		frmCodMLControls.setBounds(100, 100, 1200, 800);
		frmCodMLControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmCodMLControls.setPreferredSize(new Dimension(frmCodMLControls.getBounds().width, frmCodMLControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmCodMLControls.getPreferredSize());
		frmCodMLControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmCodMLControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][pref!,grow][pref!,grow][][]", "[][][][][][][][][][][][][][][][][][][][][][][]"));

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
		txtInputFile.setText("infile");
		panel.add(txtInputFile, "cell 1 0 4 1,growx");
		
		btnInputTree = new JButton("Input Tree");
		btnInputTree.setEnabled(false);
		btnInputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputTree);
			}
		});
		btnInputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnInputTree, "cell 0 1,growx");
		
		txtInputTree = new JTextField();
		txtInputTree.setEnabled(false);
		txtInputTree.setText("intree");
		txtInputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInputTree, "cell 1 1 4 1,growx");
			
		btnWeightFile = new JButton("Weights File");
		btnWeightFile.setEnabled(false);
		btnWeightFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightFile);
			}
		});
		btnWeightFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnWeightFile, "cell 0 2,growx");
		
		txtWeightFile = new JTextField();
		txtWeightFile.setEnabled(false);
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtWeightFile.setText("weightfile");
		panel.add(txtWeightFile, "cell 1 2 4 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputFile, "cell 0 3,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 3 4 1,growx");
		
		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputTree);
			}
		});
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputTree, "cell 0 4,growx");
		
		txtOutputTree = new JTextField();
		txtOutputTree.setText("outtree");
		txtOutputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutputTree, "cell 1 4 4 1,growx");
		
		lblSearchBest = new JLabel("Search for best tree:");
		lblSearchBest.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSearchBest.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSearchBest, "flowx,cell 0 5 2 1,alignx right");
		
		cmbxTreeSearchMethod = new JComboBox();
		cmbxTreeSearchMethod.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxTreeSearchMethod.setModel(new DefaultComboBoxModel(new String[] {"Yes", "No, use user trees in input file", "Yes, rearrange on user tree"}));
		cmbxTreeSearchMethod.setSelectedIndex(0);
		cmbxTreeSearchMethod.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(cmbxTreeSearchMethod.getSelectedIndex());
			}
		});
		panel.add(cmbxTreeSearchMethod, "cell 2 5,growx");
		
		lblUseLengths = new JLabel("Use lengths from user trees:");
		lblUseLengths.setEnabled(false);
		lblUseLengths.setHorizontalAlignment(SwingConstants.RIGHT);
		lblUseLengths.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseLengths, "flowx,cell 0 6 2 1,alignx right");
		
		rdbtnUseLengthsYes = new JRadioButton("Yes");
		rdbtnUseLengthsYes.setEnabled(false);
		rdbtnUseLengthsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseLengthsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UseLengthsToggle(true);
			}
		});
		rdbtnUseLengthsYes.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnUseLengthsYes, "cell 2 6");
		
		rdbtnUseLengthsNo = new JRadioButton("No");
		rdbtnUseLengthsNo.setSelected(true);
		rdbtnUseLengthsNo.setEnabled(false);
		rdbtnUseLengthsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseLengthsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UseLengthsToggle(false);
			}
		});
		rdbtnUseLengthsNo.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnUseLengthsNo, "cell 2 6");
		
		lblReadInputData = new JLabel("Read input data as:");
		lblReadInputData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblReadInputData, "cell 1 7,alignx right");
		
		rdbtnNucleotides = new JRadioButton("Nucleotides");
		rdbtnNucleotides.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnNucleotides.setSelected(true);
		rdbtnNucleotides.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				NucleotidesToggle(true);
			}
		});
		panel.add(rdbtnNucleotides, "cell 2 7");
		
		rdbtnAminoAcids = new JRadioButton("Amino Acids");
		rdbtnAminoAcids.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAminoAcids.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				NucleotidesToggle(false);
			}
		});
		panel.add(rdbtnAminoAcids, "cell 2 7");
		
		
		lblNucSubModel = new JLabel("Nucleotide substitiution model:");
		lblNucSubModel.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNucSubModel, "cell 1 8,alignx trailing");
		
		cmbxNucSubMethod = new JComboBox();
		cmbxNucSubMethod.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxNucSubMethod.setModel(new DefaultComboBoxModel(new String[] {"F84", "HKY", "Kimura 2-parameter", "Jukes-Cantor"}));
		cmbxNucSubMethod.setSelectedIndex(0);
		cmbxNucSubMethod.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				NucSubMethodToggle(cmbxNucSubMethod.getSelectedIndex());
			}
		});
		panel.add(cmbxNucSubMethod, "cell 2 8,growx");

		lblTTratio = new JLabel("Transition/transversion ratio:");
		lblTTratio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblTTratio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTTratio, "flowx,cell 0 9 2 1,alignx right");
		
		txtTTratio = new JTextField();
		txtTTratio.setFont(new Font("Arial", Font.PLAIN, 13));
		txtTTratio.setText("2.00");
		txtTTratio.setColumns(6);
		panel.add(txtTTratio, "cell 2 9");
		
		lblEmpBF = new JLabel("Use empirical base frequencies:");
		lblEmpBF.setHorizontalAlignment(SwingConstants.RIGHT);
		lblEmpBF.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblEmpBF, "flowx,cell 0 10 2 1,alignx right");
		
		rdbtnEmpBFYes = new JRadioButton("Yes");
		rdbtnEmpBFYes.setSelected(true);
		rdbtnEmpBFYes.setBackground(new Color(204, 255, 255));
		rdbtnEmpBFYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnEmpBFYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EmpBFToggle(true);
			}
		});
		panel.add(rdbtnEmpBFYes, "cell 2 10");
		
		rdbtnEmpBFNo = new JRadioButton("No");
		rdbtnEmpBFNo.setBackground(new Color(204, 255, 255));
		rdbtnEmpBFNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnEmpBFNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EmpBFToggle(false);
			}
		});
		panel.add(rdbtnEmpBFNo, "cell 2 10");
		
		lblBaseFreq = new JLabel("Base frequencies:");
		lblBaseFreq.setEnabled(false);
		lblBaseFreq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblBaseFreq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreq, "flowx,cell 0 11 2 1,alignx right");
		
		lblBaseFreqA = new JLabel("A:");
		lblBaseFreqA.setEnabled(false);
		lblBaseFreqA.setBackground(new Color(153, 255, 255));
		lblBaseFreqA.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqA, "cell 2 11");
		
		txtBaseFreqA = new JTextField();
		txtBaseFreqA.setEnabled(false);
		txtBaseFreqA.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqA.setText("0.25");
		txtBaseFreqA.setColumns(6);
		panel.add(txtBaseFreqA, "cell 2 11");
		
		lblBaseFreqC = new JLabel("   C:");
		lblBaseFreqC.setEnabled(false);
		lblBaseFreqC.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqC, "cell 2 11");
		
		txtBaseFreqC = new JTextField();
		txtBaseFreqC.setEnabled(false);
		txtBaseFreqC.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqC.setText("0.25");
		txtBaseFreqC.setColumns(6);
		panel.add(txtBaseFreqC, "cell 2 11");
		
		lblBaseFreqG = new JLabel("G:");
		lblBaseFreqG.setEnabled(false);
		lblBaseFreqG.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqG, "cell 2 12");
		
		txtBaseFreqG = new JTextField();
		txtBaseFreqG.setEnabled(false);
		txtBaseFreqG.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqG.setText("0.25");
		txtBaseFreqG.setColumns(6);
		panel.add(txtBaseFreqG, "cell 2 12");
		
		lblBaseFreqTU = new JLabel("T/U:");
		lblBaseFreqTU.setEnabled(false);
		lblBaseFreqTU.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqTU, "cell 2 12");
		
		txtBaseFreqTU = new JTextField();
		txtBaseFreqTU.setEnabled(false);
		txtBaseFreqTU.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqTU.setText("0.25");
		txtBaseFreqTU.setColumns(6);
		panel.add(txtBaseFreqTU, "cell 2 12");
		
		lblCodonTable = new JLabel("Genetic code to use:");
		lblCodonTable.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCodonTable, "cell 1 13,alignx trailing");
		
		cmbxCodonTable = new JComboBox();
		cmbxCodonTable.setModel(new DefaultComboBoxModel(new String[] {"Universal", "Ciliate", "Universal mitochondrial", "Vertebrate mitochondrial", "Fly mitochondrial", "Yeast mitochondrial"}));
		cmbxCodonTable.setSelectedIndex(0);
		cmbxCodonTable.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxCodonTable, "cell 2 13,growx");
		
		lblOneOmega = new JLabel("One value of omega:");
		lblOneOmega.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblOneOmega, "cell 1 14,alignx right");
		
		rdbtnOneOmegaYes = new JRadioButton("Yes");
		rdbtnOneOmegaYes.setSelected(true);
		rdbtnOneOmegaYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOneOmegaYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OneOmegaToggle(true);
			}
		});
		panel.add(rdbtnOneOmegaYes, "flowx,cell 2 14");
		
		rdbtnOneOmegaNo = new JRadioButton("No");
		rdbtnOneOmegaNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOneOmegaNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OneOmegaToggle(false);
			}
		});
		panel.add(rdbtnOneOmegaNo, "cell 2 14");

		lblNSratio = new JLabel("Nonsynonymous/synonymous ratio (omega):");
		lblNSratio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNSratio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNSratio, "flowx,cell 0 15 2 1,alignx right");
		
		txtNSratio = new JTextField();
		txtNSratio.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNSratio.setText("1.00");
		txtNSratio.setColumns(6);
		panel.add(txtNSratio, "cell 2 15");

		lblOmegasCor = new JLabel("Omegas at adjacent sites correlated:");
		lblOmegasCor.setEnabled(false);
		lblOmegasCor.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblOmegasCor, "cell 1 16,alignx right");

		rdbtnOmegasCorYes = new JRadioButton("Yes");
		rdbtnOmegasCorYes.setEnabled(false);
		rdbtnOmegasCorYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOmegasCorYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OmegaCorToggle(true);
			}
		});
		panel.add(rdbtnOmegasCorYes, "flowx,cell 2 16");
		
		rdbtnOmegasCorNo = new JRadioButton("No, independent");
		rdbtnOmegasCorNo.setEnabled(false);
		rdbtnOmegasCorNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOmegasCorNo.setSelected(true);
		rdbtnOmegasCorNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OmegaCorToggle(false);
			}
		});
		panel.add(rdbtnOmegasCorNo, "cell 2 16");
		
		lblBlockLen = new JLabel("Mean block length:");
		lblBlockLen.setEnabled(false);
		lblBlockLen.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBlockLen, "cell 2 17");
		
		txtBlockLen = new JTextField();
		txtBlockLen.setText("2");
		txtBlockLen.setEnabled(false);
		panel.add(txtBlockLen, "cell 2 17");
		txtBlockLen.setColumns(5);
		
		lblGt1 = new JLabel("(>1)");
		lblGt1.setEnabled(false);
		lblGt1.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGt1, "cell 2 17");
		
		lblCodonsWeighted = new JLabel("Codons weighted:");
		lblCodonsWeighted.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCodonsWeighted, "cell 1 18,alignx right");
		
		rdbtnCodonsWeightedYes = new JRadioButton("Yes");
		rdbtnCodonsWeightedYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCodonsWeightedYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CodonWgtToggle(true);
			}
		});
		panel.add(rdbtnCodonsWeightedYes, "flowx,cell 2 18");
		
		rdbtnCodonsWeightedNo = new JRadioButton("No");
		rdbtnCodonsWeightedNo.setSelected(true);
		rdbtnCodonsWeightedNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCodonsWeightedNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CodonWgtToggle(false);
			}
		});
		panel.add(rdbtnCodonsWeightedNo, "cell 2 18");
		
		lblOmegaValuesHMM = new JLabel("Omega values in HMM:");
		lblOmegaValuesHMM.setEnabled(false);
		lblOmegaValuesHMM.setFont(new Font("Arial", Font.BOLD, 13));
		lblOmegaValuesHMM.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblOmegaValuesHMM, "cell 1 19,alignx right");
		
		lblNumCats = new JLabel("Number of categories:");
		lblNumCats.setEnabled(false);
		lblNumCats.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNumCats, "cell 1 20,alignx trailing");
		
		cmbxNumCats = new JComboBox();
		cmbxNumCats.setEnabled(false);
		cmbxNumCats.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxNumCats.setSelectedIndex(0);
		cmbxNumCats.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxNumCats, "flowx,cell 2 20,alignx left");
		
		lblCategory = new JLabel("Category:");
		lblCategory.setEnabled(false);
		lblCategory.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCategory, "cell 2 20");
		
		cmbxCategory = new JComboBox();
		cmbxCategory.setEnabled(false);
		cmbxCategory.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxCategory.setSelectedIndex(0);
		cmbxCategory.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxCategory.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplayHMMValues(Integer.parseInt(cmbxCategory.getSelectedItem().toString()));
			}
		});
		panel.add(cmbxCategory, "cell 2 20");
		
		lblRate = new JLabel("Rate:");
		lblRate.setEnabled(false);
		lblRate.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRate, "flowx,cell 2 21");
		
		txtOmegaVal = new JTextField();
		txtOmegaVal.setEnabled(false);
		txtOmegaVal.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOmegaVal.setText("1.0");
		panel.add(txtOmegaVal, "cell 2 21");
		txtOmegaVal.setColumns(4);
		
		lblProbability = new JLabel("Probability:");
		lblProbability.setEnabled(false);
		lblProbability.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblProbability, "cell 2 21");
		
		txtOmegaProb = new JTextField();
		txtOmegaProb.setEnabled(false);
		txtOmegaProb.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOmegaProb.setText("1.0");
		panel.add(txtOmegaProb, "cell 2 21");
		txtOmegaProb.setColumns(4);

		/*** second column ***/	
		
		lblSpeedyAnalysis = new JLabel("Speedier but rougher analysis:");
		lblSpeedyAnalysis.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSpeedyAnalysis, "cell 3 5,alignx right");
		
		rdbtnSpeedyYes = new JRadioButton("Yes");
		rdbtnSpeedyYes.setSelected(true);
		rdbtnSpeedyYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSpeedyYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SpeedyToggle(true);
			}
		});
		panel.add(rdbtnSpeedyYes, "flowx,cell 4 5");
		
		rdbtnSpeedyNo = new JRadioButton("No");
		rdbtnSpeedyNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSpeedyNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SpeedyToggle(false);
			}
		});
		panel.add(rdbtnSpeedyNo, "cell 4 5");
		
		lblGlobalRearr = new JLabel("Global rearrangements:");
		lblGlobalRearr.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGlobalRearr, "cell 3 6,alignx right");
		
		rdbtnGlobalRearrYes = new JRadioButton("Yes");
		rdbtnGlobalRearrYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnGlobalRearrYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalRearrToggle(true);
			}
		});
		panel.add(rdbtnGlobalRearrYes, "flowx,cell 4 6");
		
		rdbtnGlobalRearrNo = new JRadioButton("No");
		rdbtnGlobalRearrNo.setSelected(true);
		rdbtnGlobalRearrNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnGlobalRearrNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalRearrToggle(false);
			}
		});
		panel.add(rdbtnGlobalRearrNo, "cell 4 6");

		lblRandInSeqOrd = new JLabel("Randomize input order of sequence:");
		lblRandInSeqOrd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandInSeqOrd, "cell 3 7,alignx right");
		
		rdbtnRandInSeqOrdYes = new JRadioButton("Yes");
		rdbtnRandInSeqOrdYes.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(rdbtnRandInSeqOrdYes, "flowx,cell 4 7");
		rdbtnRandInSeqOrdYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(true);
			}
		});
	
		rdbtnRandInSeqOrdNo = new JRadioButton("No, use input order");
		rdbtnRandInSeqOrdNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRandInSeqOrdNo.setSelected(true);
		panel.add(rdbtnRandInSeqOrdNo, "cell 4 7");
		rdbtnRandInSeqOrdNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(false);
			}
		});

		lblRandSeed = new JLabel("Random number seed:");
		lblRandSeed.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandSeed.setFont(new Font("Arial", Font.BOLD, 13));
		lblRandSeed.setEnabled(false);
		panel.add(lblRandSeed, "cell 3 8,alignx right");
		
		txtRandSeed = new JTextField();
		txtRandSeed.setText("1");
		txtRandSeed.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRandSeed.setEnabled(false);
		txtRandSeed.setColumns(6);
		panel.add(txtRandSeed, "cell 4 8");
		
		lblRandOdd = new JLabel("(must be odd)");
		lblRandOdd.setEnabled(false);
		lblRandOdd.setHorizontalAlignment(SwingConstants.LEFT);
		lblRandOdd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOdd, "cell 4 8");

		lblNumberJumble = new JLabel("Number of times to jumble:");
		lblNumberJumble.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberJumble.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumberJumble.setEnabled(false);
		panel.add(lblNumberJumble, "cell 3 9,alignx right");

		txtNumberJumble = new JTextField();
		txtNumberJumble.setColumns(6);
		txtNumberJumble.setText("1");
		txtNumberJumble.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberJumble.setEnabled(false);
		txtNumberJumble.setBounds(826, 248, 86, 20);
		panel.add(txtNumberJumble, "cell 4 9");
	
		lblOutRoot = new JLabel("Outgroup root:");
		lblOutRoot.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutRoot.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblOutRoot, "cell 3 10,alignx right");
		
		rdbtnOutRootYes = new JRadioButton("Yes");
		rdbtnOutRootYes.setBackground(new Color(204, 255, 255));
		rdbtnOutRootYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutRootYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutRootToggle(false);
			}
		});
		rdbtnOutRootYes.setSelected(false);
		panel.add(rdbtnOutRootYes, "cell 4 10");
		
		rdbtnOutRootNo = new JRadioButton("No, use as outgroup species");
		rdbtnOutRootNo.setBackground(new Color(204, 255, 255));
		rdbtnOutRootNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutRootNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutRootToggle(true);
			}
		});
		rdbtnOutRootNo.setSelected(true);
		panel.add(rdbtnOutRootNo, "cell 4 10");
		
		lblOutRootNum = new JLabel("Number of the outgroup:");
		lblOutRootNum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutRootNum.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutRootNum.setEnabled(false);
		panel.add(lblOutRootNum, "cell 3 11,alignx right");
		
		txtOutRootNum = new JTextField();
		txtOutRootNum.setText("1");
		txtOutRootNum.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutRootNum.setEnabled(false);
		txtOutRootNum.setColumns(6);
		panel.add(txtOutRootNum, "cell 4 11");
		
		lblAnalyzeMultData = new JLabel("Analyze multiple data sets:");
		lblAnalyzeMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAnalyzeMultData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAnalyzeMultData, "cell 3 12,alignx right");
		
		rdbtnAnalyzeMultDataYes = new JRadioButton("Yes");
		rdbtnAnalyzeMultDataYes.setBackground(new Color(204, 255, 255));
		rdbtnAnalyzeMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAnalyzeMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		rdbtnAnalyzeMultDataYes.setSelected(false);
		panel.add(rdbtnAnalyzeMultDataYes, "cell 4 12");
		
		rdbtnAnalyzeMultDataNo = new JRadioButton("No");
		rdbtnAnalyzeMultDataNo.setBackground(new Color(204, 255, 255));
		rdbtnAnalyzeMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAnalyzeMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		rdbtnAnalyzeMultDataNo.setSelected(true);
		panel.add(rdbtnAnalyzeMultDataNo, "cell 4 12");
		
		lblMultData = new JLabel("Multiple data sets or multiple weights:");
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultData.setEnabled(false);
		panel.add(lblMultData, "cell 3 13,alignx right");
		
		rdbtnDataSets = new JRadioButton("Data sets");
		rdbtnDataSets.setBackground(new Color(204, 255, 255));
		rdbtnDataSets.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(true);
			}
		});
		rdbtnDataSets.setSelected(true);
		rdbtnDataSets.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDataSets.setEnabled(false);
		panel.add(rdbtnDataSets, "cell 4 13");
		
		rdbtnWeights = new JRadioButton("Weights");
		rdbtnWeights.setBackground(new Color(204, 255, 255));
		rdbtnWeights.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(false);
			}
		});		
		rdbtnWeights.setSelected(false);
		rdbtnWeights.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWeights.setEnabled(false);
		panel.add(rdbtnWeights, "cell 4 13");

		lblHowManyData = new JLabel("Number:");
		lblHowManyData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHowManyData.setFont(new Font("Arial", Font.BOLD, 13));
		lblHowManyData.setEnabled(false);
		panel.add(lblHowManyData, "cell 4 14");
		
		txtNumSeqs = new JTextField();
		txtNumSeqs.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumSeqs.setText("1");
		txtNumSeqs.setEnabled(false);
		txtNumSeqs.setColumns(6);
		panel.add(txtNumSeqs, "cell 4 14");

		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblInputSeq, "cell 3 15,alignx right");
		
		rdbtnInputFileSeqYes = new JRadioButton("Interleaved");
		rdbtnInputFileSeqYes.setBackground(new Color(204, 255, 255));
		rdbtnInputFileSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputFileSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		rdbtnInputFileSeqYes.setSelected(true);
		panel.add(rdbtnInputFileSeqYes, "cell 4 15");
		
		rdbtnInputFileSeqNo = new JRadioButton("Sequential");
		rdbtnInputFileSeqNo.setBackground(new Color(204, 255, 255));
		rdbtnInputFileSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputFileSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		rdbtnInputFileSeqNo.setSelected(false);
		panel.add(rdbtnInputFileSeqNo, "cell 4 15");
		
		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintData, "cell 3 16,alignx right");
		
		rdbtnPrintDataYes = new JRadioButton("Yes");
		rdbtnPrintDataYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		rdbtnPrintDataYes.setSelected(false);
		panel.add(rdbtnPrintDataYes, "cell 4 16");
		
		rdbtnPrintDataNo = new JRadioButton("No");
		rdbtnPrintDataNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rdbtnPrintDataNo.setSelected(true);
		panel.add(rdbtnPrintDataNo, "cell 4 16");

		lblPrintTree = new JLabel("Print out tree:");
		lblPrintTree.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintTree, "cell 3 17,alignx right");
		
		rdbtnPrintTreeYes = new JRadioButton("Yes");
		rdbtnPrintTreeYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(true);
			}
		});
		rdbtnPrintTreeYes.setSelected(true);
		panel.add(rdbtnPrintTreeYes, "cell 4 17");
		
		rdbtnPrintTreeNo = new JRadioButton("No");
		rdbtnPrintTreeNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(false);
			}
		});
		rdbtnPrintTreeNo.setSelected(false);
		panel.add(rdbtnPrintTreeNo, "cell 4 17");
		
		lblWriteTree = new JLabel("Write out trees onto tree file:");
		lblWriteTree.setHorizontalAlignment(SwingConstants.RIGHT);
		lblWriteTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblWriteTree, "cell 3 18,alignx right");
		
		rdbtnWriteTreeYes = new JRadioButton("Yes");
		rdbtnWriteTreeYes.setBackground(new Color(204, 255, 255));
		rdbtnWriteTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(true);
			}
		});
		rdbtnWriteTreeYes.setSelected(true);
		panel.add(rdbtnWriteTreeYes, "cell 4 18");
		
		rdbtnWriteTreeNo = new JRadioButton("No");
		rdbtnWriteTreeNo.setBackground(new Color(204, 255, 255));
		rdbtnWriteTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(false);
			}
		});
		rdbtnWriteTreeNo.setSelected(false);
		panel.add(rdbtnWriteTreeNo, "cell 4 18");

		lblDotDiff = new JLabel("Use dot-differencing to display them:");
		lblDotDiff.setFont(new Font("Arial", Font.BOLD, 13));
		lblDotDiff.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblDotDiff, "cell 3 19,alignx right");

		rdbtnDotDiffYes = new JRadioButton("Yes");
		rdbtnDotDiffYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(true);
			}
		});
		rdbtnDotDiffYes.setSelected(true);
		panel.add(rdbtnDotDiffYes, "cell 4 19");

		rdbtnDotDiffNo = new JRadioButton("No");
		rdbtnDotDiffNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(false);
			}
		});
		panel.add(rdbtnDotDiffNo, "cell 4 19");
		
		lblRecHypo = new JLabel("Reconstruct hypothetical sequences:");
		lblRecHypo.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRecHypo.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRecHypo, "cell 3 20,alignx right");
		
		rdbtnRecHypoYes = new JRadioButton("Yes");
		rdbtnRecHypoYes.setBackground(new Color(204, 255, 255));
		rdbtnRecHypoYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRecHypoYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RecHypoToggle(false);
			}
		});
		rdbtnRecHypoYes.setSelected(false);
		panel.add(rdbtnRecHypoYes, "cell 4 20");
		
		rdbtnRecHypoNo = new JRadioButton("No");
		rdbtnRecHypoNo.setBackground(new Color(204, 255, 255));
		rdbtnRecHypoNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRecHypoNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RecHypoToggle(true);
			}
		});
		rdbtnRecHypoNo.setSelected(true);
		panel.add(rdbtnRecHypoNo, "cell 4 20");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInd, "cell 3 21,alignx right");
		
		rdbtnPrintIndYes = new JRadioButton("Yes");
		rdbtnPrintIndYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rdbtnPrintIndYes.setSelected(true);
		panel.add(rdbtnPrintIndYes, "cell 4 21");
		
		rdbtnPrintIndNo = new JRadioButton("No");
		rdbtnPrintIndNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rdbtnPrintIndNo.setSelected(false);
		panel.add(rdbtnPrintIndNo, "cell 4 21");
		
		btnExecute = new JButton("Execute");
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 4 22,alignx center");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				// catch last array data entries
				OmegaVal[Integer.parseInt(cmbxCategory.getSelectedItem().toString())-1]  = Double.parseDouble(txtOmegaVal.getText());
				OmegaProb[Integer.parseInt(cmbxCategory.getSelectedItem().toString())-1]  = Double.parseDouble(txtOmegaProb.getText());

				inputvals = new CodMLData();
				inputvals.infile = txtInputFile.getText();
				inputvals.intree = txtInputTree.getText();
				inputvals.wgtsfile = txtWeightFile.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.outtree = txtOutputTree.getText();
				inputvals.outtreeopt = "w";
				inputvals.TreeUseMethod = cmbxTreeSearchMethod.getSelectedItem().toString();
				inputvals.UseLengths = rdbtnUseLengthsYes.isSelected();
				inputvals.InputNuc = rdbtnNucleotides.isSelected();
				inputvals.NucSubModel = cmbxNucSubMethod.getSelectedItem().toString();
				inputvals.TTratio = Double.parseDouble(txtTTratio.getText());
				inputvals.NSratio = Double.parseDouble(txtNSratio.getText());
				inputvals.useEmpBF = rdbtnEmpBFYes.isSelected();
				inputvals.BaseFreqA = Double.parseDouble(txtBaseFreqA.getText());
				inputvals.BaseFreqC = Double.parseDouble(txtBaseFreqC.getText());
				inputvals.BaseFreqG = Double.parseDouble(txtBaseFreqG.getText());
				inputvals.BaseFreqTU = Double.parseDouble(txtBaseFreqTU.getText());				
				if (cmbxCodonTable.getSelectedIndex() == 0)
				{
					inputvals.GeneCode = "Universal";
				}
				else if (cmbxCodonTable.getSelectedIndex() == 1)
				{
					inputvals.GeneCode = "Ciliate";
				}
				else if (cmbxCodonTable.getSelectedIndex() == 2)
				{
					inputvals.GeneCode = "Mitochondrial";
				}
				else if (cmbxCodonTable.getSelectedIndex() == 3)
				{
					inputvals.GeneCode = "Vertebrate";
				}
				else if (cmbxCodonTable.getSelectedIndex() == 4)
				{
					inputvals.GeneCode = "Fly";
				}
				else // if (cmbxCodonTable.getSelectedIndex() == 5)
				{
					inputvals.GeneCode = "Yeast";
				}
				//inputvals.RateVar = cmbxRateSite.getSelectedItem().toString();
				inputvals.OneOmega = rdbtnOneOmegaYes.isSelected();
				inputvals.AdjOmegasCor = rdbtnOmegasCorYes.isSelected();
				inputvals.OmegaBlockLen = Integer.parseInt(txtBlockLen.getText());;
				inputvals.CodonsWgted = rdbtnCodonsWeightedYes.isSelected();
				inputvals.NumOmegas = Integer.parseInt(cmbxNumCats.getSelectedItem().toString());
				
				inputvals.OmegaVal1 = OmegaVal[0];
				inputvals.OmegaVal2 = OmegaVal[1];
				inputvals.OmegaVal3 = OmegaVal[2];
				inputvals.OmegaVal4 = OmegaVal[3];
				inputvals.OmegaVal5 = OmegaVal[4];
				inputvals.OmegaVal6 = OmegaVal[5];
				inputvals.OmegaVal7 = OmegaVal[6];
				inputvals.OmegaVal8 = OmegaVal[7];
				inputvals.OmegaVal9 = OmegaVal[8];
				
				inputvals.OmegaProb1 = OmegaProb[0];
				inputvals.OmegaProb2 = OmegaProb[1];
				inputvals.OmegaProb3 = OmegaProb[2];
				inputvals.OmegaProb4 = OmegaProb[3];
				inputvals.OmegaProb5 = OmegaProb[4];
				inputvals.OmegaProb6 = OmegaProb[5];
				inputvals.OmegaProb7 = OmegaProb[6];
				inputvals.OmegaProb8 = OmegaProb[7];
				inputvals.OmegaProb9 = OmegaProb[8];
				
				inputvals.SpeedAn = rdbtnSpeedyYes.isSelected();
				inputvals.GlobalRe = rdbtnGlobalRearrYes.isSelected();
				inputvals.RandInput = rdbtnRandInSeqOrdYes.isSelected();
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
				
				btnExecute.setEnabled(false);	
				String title = "Codml Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread codMLThread = new Thread() {
						public void run() {
							runCodMLThreads();
						}
			  	    };
			  	    codMLThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmCodMLControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 4 22");		
		/*** end second column ***/
		
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
		
		if (!inputvals.MultDSet){
			if (!test.FileAvailable(inputvals.wgtsfile, "Weights"))
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
		
		if (inputvals.OmegaBlockLen <= 0) {
			String msg1 = "Input value: Omega Mean Block Length must be greater than 1.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.TTratio < 0) {
			String msg1 = "Input value: Transition / Transversion Ratio cannot be negative.";
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
		double probsum = inputvals.OmegaProb1;
		
		if (i<inputvals.NumOmegas){
    		probsum += inputvals.OmegaProb2;
		}
		else
		{
			inputvals.OmegaProb2 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumOmegas){
    		probsum += inputvals.OmegaProb3;
		}
		else
		{
			inputvals.OmegaProb3 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumOmegas){
    		probsum += inputvals.OmegaProb4;
		}
		else
		{
			inputvals.OmegaProb4 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumOmegas){
    		probsum += inputvals.OmegaProb5;
		}
		else
		{
			inputvals.OmegaProb5 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumOmegas){
    		probsum += inputvals.OmegaProb6;
		}
		else
		{
			inputvals.OmegaProb6 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumOmegas){
    		probsum += inputvals.OmegaProb7;
		}
		else
		{
			inputvals.OmegaProb7 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumOmegas){
    		probsum += inputvals.OmegaProb8;
		}
		else
		{
			inputvals.OmegaProb8 = 0.0;
		}
		i++;
		
		if (i<inputvals.NumOmegas){
    		probsum += inputvals.OmegaProb9;
		}
		else
		{
			inputvals.OmegaProb9 = 0.0;
		}
		
		inputvals.OmegaProb1 = inputvals.OmegaProb1/probsum;
		inputvals.OmegaProb2 = inputvals.OmegaProb2/probsum;
		inputvals.OmegaProb3 = inputvals.OmegaProb3/probsum;
		inputvals.OmegaProb4 = inputvals.OmegaProb4/probsum;
		inputvals.OmegaProb5 = inputvals.OmegaProb5/probsum;
		inputvals.OmegaProb6 = inputvals.OmegaProb6/probsum;
		inputvals.OmegaProb7 = inputvals.OmegaProb7/probsum;
		inputvals.OmegaProb8 = inputvals.OmegaProb8/probsum;
		inputvals.OmegaProb9 = inputvals.OmegaProb9/probsum;
		
		return true;
	}
	  
	protected void runCodMLThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("codml", CodML.class);
 		}
    	catch(UnsatisfiedLinkError e)
	    {
     		new TestFileNames().LibraryMissing("CodML");
    		return;
	    }
    	
    	try
    	{
	  	    Thread codMLRunThread = new Thread() {
		  	      public void run() {
		  			// at this point we hook into the C code
		  	    	CodML codml = (CodML) Native.loadLibrary("codml", CodML.class);
		  	    	codml.codml(		
		  	        		inputvals.infile,
		  	        		inputvals.intree,
		  	        		inputvals.wgtsfile,
		  	        		inputvals.outfile,
		  	        		inputvals.outfileopt,
		  	        		inputvals.outtree,
		  	        		inputvals.outtreeopt,
		  	        		inputvals.TreeUseMethod,
		  	        		inputvals.UseLengths,
		  	        		inputvals.InputNuc,
		  	        		inputvals.NucSubModel,
		  	        		inputvals.TTratio,
		  	        		inputvals.NSratio,
		  	        		inputvals.useEmpBF,
		  	        		inputvals.BaseFreqA,
		  	        		inputvals.BaseFreqC,
		  	        		inputvals.BaseFreqG,
		  	        		inputvals.BaseFreqTU,
		  	        		inputvals.GeneCode,
		  	        		inputvals.OneOmega,
		  	        		inputvals.AdjOmegasCor,
		  	        		inputvals.OmegaBlockLen,
		  	        		inputvals.CodonsWgted,
		  	        		inputvals.NumOmegas,
		  	        		inputvals.OmegaVal1,
		  	                inputvals.OmegaVal2,
		  	                inputvals.OmegaVal3,
		  	                inputvals.OmegaVal4,
		  	                inputvals.OmegaVal5,
		  	                inputvals.OmegaVal6,
		  	                inputvals.OmegaVal7,
		  	                inputvals.OmegaVal8,
		  	                inputvals.OmegaVal9,
		  	        		inputvals.OmegaProb1,
		  	                inputvals.OmegaProb2,
		  	                inputvals.OmegaProb3,
		  	                inputvals.OmegaProb4,
		  	                inputvals.OmegaProb5,
		  	                inputvals.OmegaProb6,
		  	                inputvals.OmegaProb7,
		  	                inputvals.OmegaProb8,
		  	                inputvals.OmegaProb9,
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
	  	    codMLRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (codMLRunThread.isAlive());
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
