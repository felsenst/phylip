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

public class ProMLUserInterface {
	
	public interface ProML extends Library {
	    public void proml(
    		String infile,
    		String intree,
    		String wgtsfile,
    		String catsfile,
    		String outfile,
    		String outfileopt,
    		String outtree,
    		String outtreeopt,
    		String TreeUseMethod,
    		String ProbModel,
    		boolean UseLengths,
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
    		
    		// these are also explicitly named because JNA doesn't pass arrays gracefully
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
            
    		// these too are explicitly named because JNA doesn't pass arrays gracefully
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
    		boolean WriteTree,
    		boolean RecHypo
	        );
	    }

	public class ProMLData {
		String infile;
		String intree;
		String wgtsfile;
		String catsfile;
		String outfile;
		String outfileopt;
		String outtree;
		String outtreeopt;
		String TreeUseMethod;
		String ProbModel;
		boolean UseLengths;
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
		boolean WriteTree;
		boolean RecHypo;
	}

	private double[] SiteRateValues;
	private double[] HMMRateValues;
	private double[] HMMProbValues;
	private int lastsiteratecat;
	private int lasthmmcat;
	
	private ProMLData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;

	private JFrame frmProMLControls;
	private JButton btnInputFile;
	private JTextField txtInputFile;
	private JButton btnOutputFile;
	private JTextField txtOutputFile;
	private JLabel lblSearchBest;
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
	private JTextField txtOutputTree;
	private JButton btnCatFile;
	private JTextField txtCatFile;
	private JLabel lblRate;
	private JTextField txtSiteRate;
	private JLabel lblCatNum;
	private JComboBox cmbxRateCatNum;
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
	private JLabel lblProbMdl;
	private JComboBox cmbxProbMdl;
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
					ProMLUserInterface window = new ProMLUserInterface(args);
					window.frmProMLControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	protected void ChooseFile(JTextField file) {
		JFileChooser fileChooser = new JFileChooser(filedir);

		int option = fileChooser.showOpenDialog(frmProMLControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}

	protected void CatToggle(boolean isOneCat){
		if (isOneCat){
			rdbtnCatYes.setSelected(true);
			rdbtnCatNo.setSelected(false);
			lblNumCat.setEnabled(false);
			cmbxNumCat.setEnabled(false);
			lblCatNum.setEnabled(false);
			cmbxRateCatNum.setEnabled(false);
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
			lblCatNum.setEnabled(true);
			cmbxRateCatNum.setEnabled(true);
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
		}
		else{
			rdbtnSitesWeightYes.setSelected(true);
			rdbtnSitesWeightNo.setSelected(false);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
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
	
	protected void AnalyzeMultToggle(boolean isMult){
		if(isMult){
			rdbtnAnalyzeMultDataNo.setSelected(false);
			rdbtnAnalyzeMultDataYes.setSelected(true);
			lblMultData.setEnabled(true);
			rdbtnDataSets.setEnabled(true);
			rdbtnWeights.setEnabled(true);
			lblHowManyData.setEnabled(true);
			txtNumSeqs.setEnabled(true);
			/* 
			 * does this make sense?
			if(rdbtnDataSets.isSelected()){
				RandOrderToggle(true);
			}
			else
			{
				RandOrderToggle(false);
			}
			*/
		}
		else{
			rdbtnAnalyzeMultDataNo.setSelected(true);
			rdbtnAnalyzeMultDataYes.setSelected(false);
			lblMultData.setEnabled(false);
			rdbtnDataSets.setEnabled(false);
			rdbtnWeights.setEnabled(false);
			lblHowManyData.setEnabled(false);
			txtNumSeqs.setEnabled(false);
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
			txtOutputTree.setEnabled(true);
		}
		else{
			rdbtnWriteTreeYes.setSelected(false);
			rdbtnWriteTreeNo.setSelected(true);
			btnOutputTree.setEnabled(false);
			txtOutputTree.setEnabled(false);
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
		}
		else{
			rdbtnDataSets.setSelected(false);
			rdbtnWeights.setSelected(true);
			RandOrderToggle(false);
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

	protected void UseLengthsToggle(boolean isUseLength) {
		if (isUseLength) {
			rdbtnUseLengthsYes.setSelected(true);
			rdbtnUseLengthsNo.setSelected(false);
		} else {
			rdbtnUseLengthsYes.setSelected(false);
			rdbtnUseLengthsNo.setSelected(true);
		}
	}
	
	protected void ProbmdlToggle(int selected){
		if(selected == 0) //Jones-Taylor-Thornton
		{	
		}
		else if(selected == 1) //Heinkoff/Tiller PMB
		{
		}
		else // Dayhoff PAM
		{
		}
	}

	/**
	 * Create the application.
	 */
	public ProMLUserInterface(String[] args) {
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

		frmProMLControls = new JFrame();
		frmProMLControls.setBackground(new Color(204, 255, 255));
		frmProMLControls.setTitle("Proml");
		frmProMLControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmProMLControls.setBounds(100, 100, 1150, 750);
		frmProMLControls.setPreferredSize(new Dimension(frmProMLControls.getBounds().width, frmProMLControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmProMLControls.getPreferredSize());
		frmProMLControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmProMLControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow][pref!,grow][pref!,grow][pref!,grow]", "[][][][]"));
		
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
		panel.add(txtInputFile, "cell 1 0 5 1,growx");
		
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
		panel.add(txtInputTree, "cell 1 1 5 1,growx");
	
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
		txtWeightFile.setText("weightfile");
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtWeightFile, "cell 1 2 5 1,growx");
		
		btnCatFile = new JButton("Categories File");
		btnCatFile.setEnabled(false);
		btnCatFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtCatFile);
			}
		});
		btnCatFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnCatFile, "cell 0 3,growx");
		
		txtCatFile = new JTextField();
		txtCatFile.setEnabled(false);
		txtCatFile.setText("catfile");
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
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 4 5 1,growx");

		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputTree);
			}
		});
		panel.add(btnOutputTree, "cell 0 5,growx");
		
		txtOutputTree = new JTextField();
		txtOutputTree.setText("outtree");
		txtOutputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutputTree, "cell 1 5 5 1,growx");
		
		lblSearchBest = new JLabel("Search for best tree:");
		lblSearchBest.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSearchBest.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSearchBest, "flowx,cell 0 6 2 1,alignx right");
		
		TreeSearchMethod = new JComboBox();
		TreeSearchMethod.setModel(new DefaultComboBoxModel(new String[] {"Yes", "No, use user trees in input file", "Yes, rearrange on user tree"}));
		TreeSearchMethod.setSelectedIndex(0);
		TreeSearchMethod.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(TreeSearchMethod.getSelectedIndex());
			}
		});
		panel.add(TreeSearchMethod, "cell 2 6 2 1, growx");
		
		lblUseLengths = new JLabel("Use lengths from user trees:");
		lblUseLengths.setEnabled(false);
		lblUseLengths.setHorizontalAlignment(SwingConstants.RIGHT);
		lblUseLengths.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseLengths, "flowx,cell 0 7 2 1,alignx right");
		
		rdbtnUseLengthsYes = new JRadioButton("Yes");
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
		
		lblProbMdl = new JLabel("Probability Model:");
		lblProbMdl.setHorizontalAlignment(SwingConstants.RIGHT);
		lblProbMdl.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblProbMdl, "flowx,cell 0 8 2 1,alignx right");
		
		cmbxProbMdl = new JComboBox();
		cmbxProbMdl.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxProbMdl.setModel(new DefaultComboBoxModel(new String[] {"Jones-Taylor-Thornton", "Heinkoff/Tiller PMB", "Dayhoff PAM"}));
		cmbxProbMdl.setSelectedIndex(0);
		cmbxProbMdl.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ProbmdlToggle(cmbxProbMdl.getSelectedIndex());
			}
		});
		panel.add(cmbxProbMdl, "cell 2 8 2 1, growx");

		lblCatSites = new JLabel("One category of sites:");
		lblCatSites.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCatSites.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCatSites, "flowx,cell 0 9 2 1,alignx right");
		
		rdbtnCatYes = new JRadioButton("Yes");
		rdbtnCatYes.setSelected(true);
		rdbtnCatYes.setBackground(new Color(204, 255, 255));
		rdbtnCatYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCatYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CatToggle(true);
			}
		});
		panel.add(rdbtnCatYes, "cell 2 9");
		
		rdbtnCatNo = new JRadioButton("No");
		rdbtnCatNo.setBackground(new Color(204, 255, 255));
		rdbtnCatNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCatNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CatToggle(false);
			}
		});
		panel.add(rdbtnCatNo, "cell 2 9");
		
		lblNumCat = new JLabel("Number of site categories:");
		lblNumCat.setEnabled(false);
		lblNumCat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumCat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNumCat, "flowx,cell 0 10 2 1,alignx right");
		
		cmbxNumCat = new JComboBox();
		cmbxNumCat.setEnabled(false);
		cmbxNumCat.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxNumCat.setSelectedIndex(0);
		cmbxNumCat.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxNumCat, "cell 2 10");

		lblCatNum = new JLabel("Category:");
		lblCatNum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCatNum.setFont(new Font("Arial", Font.BOLD, 13));
		lblCatNum.setEnabled(false);
		panel.add(lblCatNum, "cell 2 10,alignx right");
		
		cmbxRateCatNum = new JComboBox();
		cmbxRateCatNum.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxRateCatNum.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxRateCatNum.setEnabled(false);
		cmbxRateCatNum.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplaySiteRateValue(Integer.parseInt(cmbxRateCatNum.getSelectedItem().toString()));
			}
		});
		panel.add(cmbxRateCatNum, "cell 3 10");
		
		lblRate = new JLabel("Rate:");
		lblRate.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRate.setFont(new Font("Arial", Font.BOLD, 13));
		lblRate.setEnabled(false);
		panel.add(lblRate, "cell 3 10");
		
		txtSiteRate = new JTextField();
		txtSiteRate.setText("1.0");
		txtSiteRate.setHorizontalAlignment(SwingConstants.CENTER);
		txtSiteRate.setFont(new Font("Arial", Font.PLAIN, 13));
		txtSiteRate.setEnabled(false);
		txtSiteRate.setColumns(5);
		panel.add(txtSiteRate, "cell 3 10");
	
		lblRateSite = new JLabel("Rate variation among sites:");
		lblRateSite.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRateSite.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRateSite, "flowx,cell 0 11 2 1,alignx right");
		
		cmbxRateSite = new JComboBox();
		cmbxRateSite.setModel(new DefaultComboBoxModel(new String[] {"Constant rate", "Gamma distance rates", "Gamma + invariant sites","User-defined HMM of rates"}));
		cmbxRateSite.setSelectedIndex(0);
		cmbxRateSite.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxRateSite.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				RatesiteToggle(cmbxRateSite.getSelectedIndex());
			}
		});
		panel.add(cmbxRateSite, "cell 2 11 2 1, growx");
		
		lblRatesAdjCor = new JLabel("Rates on adjacent sites correlated:");
		lblRatesAdjCor.setEnabled(false);
		lblRatesAdjCor.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRatesAdjCor.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRatesAdjCor, "flowx,cell 0 12 2 1,alignx right");
		
		rdbtnAdjCorYes = new JRadioButton("Yes");
		rdbtnAdjCorYes.setEnabled(false);
		rdbtnAdjCorYes.setSelected(false);
		rdbtnAdjCorYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAdjCorYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AdjCorToggle(true);
			}
		});
		rdbtnAdjCorYes.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnAdjCorYes, "cell 2 12");
		
		rdbtnAdjCorNo = new JRadioButton("No");
		rdbtnAdjCorNo.setSelected(true);
		rdbtnAdjCorNo.setEnabled(false);
		rdbtnAdjCorNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAdjCorNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AdjCorToggle(false);
			}
		});
		rdbtnAdjCorNo.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnAdjCorNo, "cell 2 12");
		
		lblMeanLen = new JLabel("Mean block length:");
		lblMeanLen.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMeanLen.setFont(new Font("Arial", Font.BOLD, 13));
		lblMeanLen.setEnabled(false);
		panel.add(lblMeanLen, "cell 3 12");
		
		txtBlockLen = new JTextField();
		txtBlockLen.setText("1");
		txtBlockLen.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBlockLen.setEnabled(false);
		txtBlockLen.setColumns(5);
		panel.add(txtBlockLen, "cell 3 12");
		
		lblCoeffVar = new JLabel("Coefficient of variation:");
		lblCoeffVar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCoeffVar.setFont(new Font("Arial", Font.BOLD, 13));
		lblCoeffVar.setEnabled(false);
		panel.add(lblCoeffVar, "flowx,cell 0 13 2 1,alignx right");

		txtCoeffVar = new JTextField();
		txtCoeffVar.setEnabled(false);
		txtCoeffVar.setText("1.0");
		txtCoeffVar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCoeffVar.setColumns(6);
		panel.add(txtCoeffVar, "cell 2 13");
		
		lblCoeffVarNote = new JLabel("(for gamma dist = 1/\u221Aalpha)");
		lblCoeffVarNote.setHorizontalAlignment(SwingConstants.LEFT);
		lblCoeffVarNote.setFont(new Font("Arial", Font.BOLD, 13));
		lblCoeffVarNote.setVisible(false);
		panel.add(lblCoeffVarNote, "cell 3 13");
		
		lblHMM = new JLabel("Rates for HMM:");
		lblHMM.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMM.setFont(new Font("Arial", Font.BOLD, 13));
		lblHMM.setEnabled(false);
		panel.add(lblHMM, "flowx,cell 0 14 2 1,alignx right");
		
		lblInvar = new JLabel("");
		lblInvar.setHorizontalAlignment(SwingConstants.LEFT);
		lblInvar.setFont(new Font("Arial", Font.BOLD, 13));
		lblInvar.setVisible(false);
		panel.add(lblInvar, "cell 2 14");
		
		lblNumHMMcats = new JLabel("Number of rate categories:");
		lblNumHMMcats.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumHMMcats.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumHMMcats.setEnabled(false);
		panel.add(lblNumHMMcats, "flowx,cell 0 15 2 1,alignx right");

		cmbxHMMcount = new JComboBox();
		cmbxHMMcount.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxHMMcount.setSelectedIndex(0);
		cmbxHMMcount.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxHMMcount.setEnabled(false);
		panel.add(cmbxHMMcount, "cell 2 15");
		
		lblHMMCat = new JLabel("Category:");
		lblHMMCat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMMCat.setFont(new Font("Arial", Font.BOLD, 13));
		lblHMMCat.setEnabled(false);
		panel.add(lblHMMCat, "cell 3 15");
		
		cmbxHMMCat = new JComboBox();
		cmbxHMMCat.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxHMMCat.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxHMMCat.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplayHMMValues(Integer.parseInt(cmbxHMMCat.getSelectedItem().toString()));
			}
		});
		cmbxHMMCat.setEnabled(false);
		panel.add(cmbxHMMCat, "cell 3 15");
		
		lblHMMRate = new JLabel("Rate:");
		lblHMMRate.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHMMRate.setFont(new Font("Arial", Font.BOLD, 13));
		lblHMMRate.setEnabled(false);
		panel.add(lblHMMRate, "cell 2 16");
		
		txtHMMRate = new JTextField();
		txtHMMRate.setText("1.0");
		txtHMMRate.setHorizontalAlignment(SwingConstants.CENTER);
		txtHMMRate.setFont(new Font("Arial", Font.PLAIN, 13));
		txtHMMRate.setEnabled(false);
		txtHMMRate.setColumns(6);
		panel.add(txtHMMRate, "cell 2 16");
		
		lblProb = new JLabel("Probablilty:");
		lblProb.setHorizontalAlignment(SwingConstants.RIGHT);
		lblProb.setFont(new Font("Arial", Font.BOLD, 13));
		lblProb.setEnabled(false);
		panel.add(lblProb, "cell 3 16");
		
		txtHMMProb = new JTextField();
		txtHMMProb.setText("1.0");
		txtHMMProb.setHorizontalAlignment(SwingConstants.CENTER);
		txtHMMProb.setFont(new Font("Arial", Font.PLAIN, 13));
		txtHMMProb.setEnabled(false);
		txtHMMProb.setColumns(6);
		panel.add(txtHMMProb, "cell 3 16");
	
		lblFracInvar = new JLabel("Fraction of sites invariant:");
		lblFracInvar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblFracInvar.setFont(new Font("Arial", Font.BOLD, 13));
		lblFracInvar.setEnabled(false);
		panel.add(lblFracInvar, "flowx,cell 0 17 2 1,alignx right");
	
		txtFracInvar = new JTextField();
		txtFracInvar.setText("1.0");
		txtFracInvar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtFracInvar.setEnabled(false);
		txtFracInvar.setColumns(6);
		panel.add(txtFracInvar, "cell 2 17");
	
		lblSitesWeight = new JLabel("Sites weighted:");
		lblSitesWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSitesWeight.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSitesWeight, "flowx,cell 0 18 2 1,alignx right");
		
		rdbtnSitesWeightYes = new JRadioButton("Yes");
		rdbtnSitesWeightYes.setBackground(new Color(204, 255, 255));
		rdbtnSitesWeightYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesWeightYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SitesWeightToggle(false);
			}
		});
		rdbtnSitesWeightYes.setSelected(false);
		panel.add(rdbtnSitesWeightYes, "cell 2 18");
		
		rdbtnSitesWeightNo = new JRadioButton("No");
		rdbtnSitesWeightNo.setBackground(new Color(204, 255, 255));
		rdbtnSitesWeightNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesWeightNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SitesWeightToggle(true);
			}
		});
		rdbtnSitesWeightNo.setSelected(true);
		panel.add(rdbtnSitesWeightNo, "cell 2 18");
		
		lblSpeedAn = new JLabel("Speedier but rougher analysis:");
		lblSpeedAn.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSpeedAn.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSpeedAn, "flowx,cell 0 19 2 1,alignx right");
		
		rdbtnSpeedAnYes = new JRadioButton("Yes");
		rdbtnSpeedAnYes.setBackground(new Color(204, 255, 255));
		rdbtnSpeedAnYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSpeedAnYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SpeedAnToggle(true);
			}
		});
		rdbtnSpeedAnYes.setSelected(true);
		panel.add(rdbtnSpeedAnYes, "cell 2 19");
		
		rdbtnSpeedAnNo = new JRadioButton("No");
		rdbtnSpeedAnNo.setBackground(new Color(204, 255, 255));
		rdbtnSpeedAnNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSpeedAnNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SpeedAnToggle(false);
			}
		});
		rdbtnSpeedAnNo.setSelected(false);
		panel.add(rdbtnSpeedAnNo, "cell 2 19");

		
		/**** start of column 2 ****/
		
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
		rdbtnGlobalReYes.setSelected(false);
		panel.add(rdbtnGlobalReYes, "cell 5 6");
		
		rdbtnGlobalReNo = new JRadioButton("No");
		rdbtnGlobalReNo.setBackground(new Color(204, 255, 255));
		rdbtnGlobalReNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnGlobalReNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalReToggle(true);
			}
		});
		rdbtnGlobalReNo.setSelected(true);
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
		rdbtnRandOrderYes.setSelected(false);
		panel.add(rdbtnRandOrderYes, "cell 5 7");
		
		rdbtnRandOrderNo = new JRadioButton("No, use input order");
		rdbtnRandOrderNo.setBackground(new Color(204, 255, 255));
		rdbtnRandOrderNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRandOrderNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(false);
			}
		});
		rdbtnRandOrderNo.setSelected(true);
		panel.add(rdbtnRandOrderNo, "cell 5 7");
		
		lblRandSeed = new JLabel("Random number seed:");
		lblRandSeed.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandSeed.setFont(new Font("Arial", Font.BOLD, 13));
		lblRandSeed.setEnabled(false);
		panel.add(lblRandSeed, "cell 4 8,alignx right");
		
		txtRandSeed = new JTextField();
		txtRandSeed.setText("1");
		txtRandSeed.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRandSeed.setEnabled(false);
		txtRandSeed.setColumns(6);
		panel.add(txtRandSeed, "cell 5 8");
		
		lblRandOdd = new JLabel("(must be odd)");
		lblRandOdd.setEnabled(false);
		lblRandOdd.setHorizontalAlignment(SwingConstants.LEFT);
		lblRandOdd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOdd, "cell 5 8");

		lblNumberJumble = new JLabel("Number of times to jumble:");
		lblNumberJumble.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberJumble.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumberJumble.setEnabled(false);
		panel.add(lblNumberJumble, "cell 4 9,alignx right");

		txtNumberJumble = new JTextField();
		txtNumberJumble.setText("1");
		txtNumberJumble.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberJumble.setEnabled(false);
		txtNumberJumble.setColumns(6);
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
		rdbtnOutRootYes.setSelected(false);
		panel.add(rdbtnOutRootYes, "cell 5 10");
		
		rdbtnOutRootNo = new JRadioButton("No, use as outgroup species");
		rdbtnOutRootNo.setBackground(new Color(204, 255, 255));
		rdbtnOutRootNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutRootNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutRootToggle(true);
			}
		});
		rdbtnOutRootNo.setSelected(true);
		panel.add(rdbtnOutRootNo, "cell 5 10");
		
		lblOutRootNum = new JLabel("Number of the outgroup:");
		lblOutRootNum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutRootNum.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutRootNum.setEnabled(false);
		panel.add(lblOutRootNum, "cell 4 11,alignx right");
		
		txtOutRootNum = new JTextField();
		txtOutRootNum.setText("1");
		txtOutRootNum.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutRootNum.setEnabled(false);
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
				AnalyzeMultToggle(true);
			}
		});
		rdbtnAnalyzeMultDataYes.setSelected(false);
		panel.add(rdbtnAnalyzeMultDataYes, "cell 5 12");
		
		rdbtnAnalyzeMultDataNo = new JRadioButton("No");
		rdbtnAnalyzeMultDataNo.setBackground(new Color(204, 255, 255));
		rdbtnAnalyzeMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAnalyzeMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AnalyzeMultToggle(false);
			}
		});
		rdbtnAnalyzeMultDataNo.setSelected(true);
		panel.add(rdbtnAnalyzeMultDataNo, "cell 5 12");
		
		lblMultData = new JLabel("Multiple data sets or multiple weights:");
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultData.setEnabled(false);
		panel.add(lblMultData, "cell 4 13,alignx right");
		
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
		panel.add(rdbtnDataSets, "cell 5 13");
		
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
		panel.add(rdbtnWeights, "cell 5 13");

		lblHowManyData = new JLabel("Number:");
		lblHowManyData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHowManyData.setFont(new Font("Arial", Font.BOLD, 13));
		lblHowManyData.setEnabled(false);
		panel.add(lblHowManyData, "cell 5 14");
		
		txtNumSeqs = new JTextField();
		txtNumSeqs.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumSeqs.setText("1");
		txtNumSeqs.setEnabled(false);
		txtNumSeqs.setColumns(6);
		panel.add(txtNumSeqs, "cell 5 14");

		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblInputSeq, "cell 4 15,alignx right");
		
		rdbtnInputSeqYes = new JRadioButton("Interleaved");
		rdbtnInputSeqYes.setBackground(new Color(204, 255, 255));
		rdbtnInputSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		rdbtnInputSeqYes.setSelected(true);
		panel.add(rdbtnInputSeqYes, "cell 5 15");
		
		rdbtnInputSeqNo = new JRadioButton("Sequential");
		rdbtnInputSeqNo.setBackground(new Color(204, 255, 255));
		rdbtnInputSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		rdbtnInputSeqNo.setSelected(false);
		panel.add(rdbtnInputSeqNo, "cell 5 15");
		
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
		rdbtnPrintDataYes.setSelected(false);
		panel.add(rdbtnPrintDataYes, "cell 5 16");
		
		rdbtnPrintDataNo = new JRadioButton("No");
		rdbtnPrintDataNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rdbtnPrintDataNo.setSelected(true);
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
		rdbtnPrintTreeYes.setSelected(true);
		panel.add(rdbtnPrintTreeYes, "cell 5 17");
		
		rdbtnPrintTreeNo = new JRadioButton("No");
		rdbtnPrintTreeNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(false);
			}
		});
		rdbtnPrintTreeNo.setSelected(false);
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
		rdbtnWriteTreeYes.setSelected(true);
		panel.add(rdbtnWriteTreeYes, "cell 5 18");
		
		rdbtnWriteTreeNo = new JRadioButton("No");
		rdbtnWriteTreeNo.setBackground(new Color(204, 255, 255));
		rdbtnWriteTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(false);
			}
		});
		rdbtnWriteTreeNo.setSelected(false);
		panel.add(rdbtnWriteTreeNo, "cell 5 18");
		
		lblRecHypo = new JLabel("Reconstruct hypothetical sequences:");
		lblRecHypo.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRecHypo.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRecHypo, "cell 4 19,alignx right");
		
		rdbtnRecHypoYes = new JRadioButton("Yes");
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
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInd, "cell 4 20,alignx right");
		
		rdbtnPrintIndYes = new JRadioButton("Yes");
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
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				// catch last array data entries
				SiteRateValues[Integer.parseInt(cmbxRateCatNum.getSelectedItem().toString())-1] = Double.parseDouble(txtSiteRate.getText());
				HMMRateValues[Integer.parseInt(cmbxHMMCat.getSelectedItem().toString())-1]  = Double.parseDouble(txtHMMRate.getText());
				HMMProbValues[Integer.parseInt(cmbxHMMCat.getSelectedItem().toString())-1]  = Double.parseDouble(txtHMMProb.getText());

				inputvals = new ProMLData();
				inputvals.infile = txtInputFile.getText();
				inputvals.intree = txtInputTree.getText();
				inputvals.wgtsfile = txtWeightFile.getText();
				inputvals.catsfile = txtCatFile.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.outtree = txtOutputTree.getText();
				inputvals.outtreeopt = "w";
				inputvals.TreeUseMethod = TreeSearchMethod.getSelectedItem().toString();
				if (cmbxProbMdl.getSelectedIndex() == 0)
				{
					inputvals.ProbModel = "JTT";
				}
				else if (cmbxProbMdl.getSelectedIndex() == 1)
				{
					inputvals.ProbModel = "PMB";
				}
				else
				{
					inputvals.ProbModel = "PAM";				
				}
				inputvals.UseLengths = rdbtnUseLengthsYes.isSelected();
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
				inputvals.InputSeq = rdbtnInputSeqYes.isSelected();
				inputvals.PrintData = rdbtnPrintDataYes.isSelected();
				inputvals.PrintInd = rdbtnPrintIndYes.isSelected();
				inputvals.PrintTree = rdbtnPrintTreeYes.isSelected();
				inputvals.WriteTree = rdbtnWriteTreeYes.isSelected();
				inputvals.RecHypo = rdbtnRecHypoYes.isSelected();
				
				btnExecute.setEnabled(false);	
				String title = "Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread promlThread = new Thread() {
						public void run() {
							runPromlThreads();
						}
			  	    };
			  	    promlThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 5 21");
		
		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				if(phylipCall)
				{
					frmProMLControls.dispose();
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
	
	protected boolean checkInputVals(){
		
		// check files	
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

		// check data
		if (inputvals.BlockLen <= 0) {
			String msg1 = "Input value: Mean Block Length must be greater than 0.";
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
	
	protected void runPromlThreads() {
    	try
    	{
    		// see if library exists
      		Native.loadLibrary("proml", ProML.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("ProML");
    		return;
    	}
    	
		try 
		{
	  	    Thread proMLRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			ProML proml = (ProML) Native.loadLibrary("proml", ProML.class);		  			
		  	        proml.proml(		
		          		inputvals.infile,
		          		inputvals.intree,
		          		inputvals.wgtsfile,
		          		inputvals.catsfile,
		          		inputvals.outfile,
		          		inputvals.outfileopt,
		          		inputvals.outtree,
		          		inputvals.outtreeopt,
		          		inputvals.TreeUseMethod,
		          		inputvals.ProbModel,
		          		inputvals.UseLengths,
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
		          		inputvals.RecHypo); 
		  	    };
	  	    };
	  	  	proMLRunThread.start();

	  	  	if (inputvals.PrintInd)
	  	  	{
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (proMLRunThread.isAlive());
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