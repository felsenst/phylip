package phylip;
import java.awt.EventQueue;
import javax.swing.JFileChooser;
import java.io.File;
import javax.swing.JFrame;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import javax.swing.JLabel;
import javax.swing.JRadioButton;
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

public class ContMLUserInterface {
	 public interface ContML extends Library {
		 public void contml(
			String infile,
			String intree,
			String outfile,
			String outfileopt,
			String outtree,
			String outtreeopt,
			boolean BestTree,
			boolean UseLengths,
			boolean GeneFreq,
			boolean AllAlleles,
			boolean OutRoot,
			int OutNum,
			boolean GlobalRearr,
			boolean RandInput,
			int RandNum,
			int Njumble,
			boolean MultData,
			int NumSets,
			boolean PrintData,
			boolean PrintInd,
			boolean PrintTree,
			boolean WriteTree);
	 }
		
	public class ContMLData {
		String infile;
		String intree;
		String outfile;
		String outfileopt;
		String outtree;
		String outtreeopt;
		boolean BestTree;
		boolean UseLengths;
		boolean GeneFreq;
		boolean AllAlleles;
		boolean OutRoot;
		int OutNum;
		boolean GlobalRearr;
		boolean RandInput;
		int RandNum;
		int Njumble;
		boolean MultData;
		int NumSets;
		boolean PrintData;
		boolean PrintInd;
		boolean PrintTree;
		boolean WriteTree;
	}

	private ContMLData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;
	
	private JFrame frmContMLControls;
	private JButton btnInput;
	private JTextField txtInFile;
	private JTextField txtInTree;
	private JButton btnInputTree;
	private JButton btnOutput;
	private JTextField txtOutFile;
	private JButton btnOutputTree;
	private JTextField txtOutTree;
	private JLabel lblSearchBest;
	private JRadioButton rbtnSearchBestYes;
	private JRadioButton rbtnSearchBestNo;
	private JLabel lblUseLengths;
	private JRadioButton rdbtnUseLengthsYes;
	private JRadioButton rdbtnUseLengthsNo;
	private JRadioButton rbtnGeneFreq;
	private JRadioButton rbtnContChar;
	private JLabel lblDataType;
	private JLabel lblAlleles;
	private JRadioButton rbtnAllelYes;
	private JRadioButton rbtnAllelNo;
	private JLabel lblOutGroup;
	private JRadioButton rbtnOutGrpYes;
	private JRadioButton rbtnOutGrpNo;
	private JLabel lblOutRootNum;
	private JTextField txtOutRootNum;
	private JLabel lblGlobal;
	private JRadioButton rbtnGlobalYes;
	private JRadioButton rbtnGlobalNo;
	private JLabel lblMultSets;
	private JRadioButton rbtnMultSetsYes;
	private JRadioButton rbtnMultSetsNo;
	private JLabel lblNumSet;
	private JTextField txtNumSets;
	private JLabel lblRandOrder;
	private JRadioButton rbtnRandOrderYes;
	private JRadioButton rbtnRandOrderNo;
	private JLabel lblRandSeed;
	private JTextField txtRandSeed;
	private JLabel lblRandOdd;
	private JLabel lblNumberJumble;
	private JTextField txtNumberJumble;
	private JLabel lblPrintInput;
	private JRadioButton rbtnPrintInputYes;
	private JRadioButton rbtnPrintInputNo;
	private JLabel lblPrintTree;
	private JRadioButton rbtnPrintTreeYes;
	private JRadioButton rbtnPrintTreeNo;
	private JLabel lblTreeFile;
	private JRadioButton rbtnTreeFileYes;
	private JRadioButton rbtnTreeFileNo;
	private JLabel lblPrintInd;
	private JRadioButton rbtnPrintIndYes;
	private JRadioButton rbtnPrintIndNo;
	private JButton btnExecute;
	private JButton btnQuit;
	
	private JScrollPane scrollPane;
	private JPanel panel;
	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					ContMLUserInterface window = new ContMLUserInterface(args);
					window.frmContMLControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmContMLControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}
	
	protected void DataTypeToggle(boolean isGene){
		if(isGene){
			rbtnGeneFreq.setSelected(true);
			rbtnContChar.setSelected(false);
		}
		else
		{
			rbtnGeneFreq.setSelected(false);
			rbtnContChar.setSelected(true);
		}
	}

	protected void OutRootToggle(boolean isOut){
		if(isOut){
			rbtnOutGrpNo.setSelected(true);
			rbtnOutGrpYes.setSelected(false);
			lblOutRootNum.setEnabled(false);
			txtOutRootNum.setEnabled(false);
		}
		else{
			rbtnOutGrpNo.setSelected(false);
			rbtnOutGrpYes.setSelected(true);
			lblOutRootNum.setEnabled(true);
			txtOutRootNum.setEnabled(true);
		}
	}
	
	
	protected void AnalyzeMultToggle(boolean isMult){
		if(isMult){
			rbtnMultSetsNo.setSelected(false);
			rbtnMultSetsYes.setSelected(true);
			lblNumSet.setEnabled(true);
			txtNumSets.setEnabled(true);
		}
		else{
			rbtnMultSetsNo.setSelected(true);
			rbtnMultSetsYes.setSelected(false);
			lblNumSet.setEnabled(false);
			txtNumSets.setEnabled(false);
		}
	}
	
	protected void PrintDataToggle(boolean isData){
		if(isData){
			rbtnPrintInputYes.setSelected(true);
			rbtnPrintInputNo.setSelected(false);
		}
		else
		{
			rbtnPrintInputYes.setSelected(false);
			rbtnPrintInputNo.setSelected(true);
		}
	}
	
	protected void PrintIndToggle(boolean isInd){
		if(isInd){
			rbtnPrintIndYes.setSelected(true);
			rbtnPrintIndNo.setSelected(false);
		}
		else{
			rbtnPrintIndYes.setSelected(false);
			rbtnPrintIndNo.setSelected(true);
		}
	}
	
	protected void PrintTreeToggle(boolean isTree){
		if(isTree){
			rbtnPrintTreeYes.setSelected(true);
			rbtnPrintTreeNo.setSelected(false);
		}
		else{
			rbtnPrintTreeYes.setSelected(false);
			rbtnPrintTreeNo.setSelected(true);
		}
	}
	
	protected void WriteTreeToggle(boolean isWrite){
		if(isWrite){
			rbtnTreeFileYes.setSelected(true);
			rbtnTreeFileNo.setSelected(false);
			btnOutputTree.setEnabled(true);
			txtOutTree.setEnabled(true);
		}
		else{
			rbtnTreeFileYes.setSelected(false);
			rbtnTreeFileNo.setSelected(true);
			btnOutputTree.setEnabled(false);
			txtOutTree.setEnabled(false);
		}
	}
	
	protected void RandOrderToggle(boolean isRand){
		if(isRand){
			rbtnRandOrderNo.setSelected(false);
			rbtnRandOrderYes.setSelected(true);
			lblRandSeed.setEnabled(true);
			txtRandSeed.setEnabled(true);
			lblNumberJumble.setEnabled(true);
			txtNumberJumble.setEnabled(true);
			lblRandOdd.setEnabled(true);
		}
		else{
			rbtnRandOrderNo.setSelected(true);
			rbtnRandOrderYes.setSelected(false);
			lblRandSeed.setEnabled(false);
			txtRandSeed.setEnabled(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEnabled(false);
			lblRandOdd.setEnabled(false);
		}
	}

	protected void GlobalRearrToggle(boolean isAllowed) {
		if (isAllowed) {
			rbtnGlobalYes.setSelected(true);
			rbtnGlobalNo.setSelected(false);
		} else {
			rbtnGlobalYes.setSelected(false);
			rbtnGlobalNo.setSelected(true);
		}
	}
	
	protected void IntreeToggle(boolean useIntree){
		if(useIntree){
			txtInTree.setEnabled(true);
			btnInputTree.setEnabled(true);
			rbtnSearchBestYes.setSelected(false);
			rbtnSearchBestNo.setSelected(true);
			lblUseLengths.setEnabled(true);
			rdbtnUseLengthsYes.setEnabled(true);
			rdbtnUseLengthsNo.setEnabled(true);
		}
		else{
			txtInTree.setEnabled(false);
			btnInputTree.setEnabled(false);
			rbtnSearchBestYes.setSelected(true);
			rbtnSearchBestNo.setSelected(false);
			lblUseLengths.setEnabled(false);
			rdbtnUseLengthsYes.setEnabled(false);
			rdbtnUseLengthsNo.setEnabled(false);
		}
	}
	
	protected void AllelesToggle(boolean isYes){
		if (isYes)
		{
			rbtnAllelYes.setSelected(true);
			rbtnAllelNo.setSelected(false);
		} else {
			rbtnAllelYes.setSelected(false);
			rbtnAllelNo.setSelected(true);
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
	public ContMLUserInterface(String[] args) {
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
		
		filedir = System.getProperty("user.dir");

		frmContMLControls = new JFrame();
		frmContMLControls.setBackground(new Color(204, 255, 255));
		frmContMLControls.setFont(new Font("Arial", Font.BOLD, 13));
		frmContMLControls.setTitle("Contml");
		frmContMLControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmContMLControls.setBounds(100, 100, 680, 685);
		frmContMLControls.setPreferredSize(new Dimension(frmContMLControls.getBounds().width, frmContMLControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmContMLControls.getPreferredSize());
		frmContMLControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmContMLControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow]", "[][][][]"));
		
		btnInput = new JButton("Input File");
		btnInput.setFont(new Font("Arial", Font.BOLD, 13));
		btnInput.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInFile);
			}
		});
		panel.add(btnInput, "cell 0 0,growx");
		
		txtInFile = new JTextField();
		txtInFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInFile.setText("infile");
		txtInFile.setBounds(168, 11, 412, 20);
		panel.add(txtInFile, "cell 1 0 2 1,growx");
		
		btnInputTree = new JButton("Input Tree");
		btnInputTree.setEnabled(false);
		btnInputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInTree);
			}
		});
		btnInputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnInputTree, "cell 0 1,growx");
		
		txtInTree = new JTextField();
		txtInTree.setEnabled(false);
		txtInTree.setText("intree");
		txtInTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInTree, "cell 1 1 2 1,growx");
		
		btnOutput = new JButton("Output File");
		btnOutput.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutput.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutFile);
			}
		});
		panel.add(btnOutput, "cell 0 2,growx");
		
		txtOutFile = new JTextField();
		txtOutFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutFile.setText("outfile");
		txtOutFile.setBounds(168, 69, 412, 20);
		panel.add(txtOutFile, "cell 1 2 2 1,growx");
		

		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutTree);
			}
		});
		panel.add(btnOutputTree, "cell 0 3,growx");
		
		txtOutTree = new JTextField();
		txtOutTree.setText("outtree");
		txtOutTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutTree, "cell 1 3 2 1,growx");
		
		lblSearchBest = new JLabel("Search for best tree:");
		lblSearchBest.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSearchBest, "flowx,cell 0 4 2 1,alignx right");

		rbtnSearchBestYes = new JRadioButton("Yes");
		rbtnSearchBestYes.setSelected(true);
		rbtnSearchBestYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSearchBestYes.setBackground(new Color(204, 255, 255));
		rbtnSearchBestYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(false);
			}
		});
		panel.add(rbtnSearchBestYes, "cell 2 4");
		
		rbtnSearchBestNo = new JRadioButton("No, use user input tree");
		rbtnSearchBestNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSearchBestNo.setBackground(new Color(204, 255, 255));
		rbtnSearchBestNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(true);
			}
		});
		panel.add(rbtnSearchBestNo, "cell 2 4 2 1");	
		
		lblUseLengths = new JLabel("Use lengths from user trees:");
		lblUseLengths.setEnabled(false);
		lblUseLengths.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseLengths, "flowx,cell 0 5 2 1,alignx right");
		
		rdbtnUseLengthsYes = new JRadioButton("Yes");
		rdbtnUseLengthsYes.setEnabled(false);
		rdbtnUseLengthsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseLengthsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UseLengthsToggle(true);
			}
		});
		rdbtnUseLengthsYes.setBackground(new Color(204, 255, 255));
		panel.add(rdbtnUseLengthsYes, "cell 2 5");
		
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
		panel.add(rdbtnUseLengthsNo, "cell 2 5");

		
		lblDataType = new JLabel("Input data type:");
		lblDataType.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDataType, "flowx,cell 0 6 2 1,alignx right");

		rbtnGeneFreq = new JRadioButton("Gene frequencies");
		rbtnGeneFreq.setSelected(true);
		rbtnGeneFreq.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnGeneFreq.setBackground(new Color(204, 255, 255));
		rbtnGeneFreq.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DataTypeToggle(true);
			}
		});
		panel.add(rbtnGeneFreq, "cell 2 6");
		
		rbtnContChar = new JRadioButton("Continuous characters");
		rbtnContChar.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnContChar.setBackground(new Color(204, 255, 255));
		rbtnContChar.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DataTypeToggle(false);
			}
		});
		panel.add(rbtnContChar, "cell 2 6");
		
		lblAlleles = new JLabel("All alleles present at each locus:");
		lblAlleles.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAlleles, "flowx,cell 0 7 2 1,alignx right");
		
		rbtnAllelYes = new JRadioButton("Yes");
		rbtnAllelYes.setSelected(false);
		rbtnAllelYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAllelYes.setBackground(new Color(204, 255, 255));
		rbtnAllelYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AllelesToggle(true);
			}
		});
		panel.add(rbtnAllelYes, "cell 2 7");
		
		rbtnAllelNo = new JRadioButton("No, one absent at each locus");
		rbtnAllelNo.setSelected(true);
		rbtnAllelNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAllelNo.setBackground(new Color(204, 255, 255));
		rbtnAllelNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AllelesToggle(false);
			}
		});
		panel.add(rbtnAllelNo, "cell 2 7 2 1");
		
		lblOutGroup = new JLabel("Outgroup root:");
		lblOutGroup.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblOutGroup, "flowx,cell 0 8 2 1,alignx right");
		
		rbtnOutGrpYes = new JRadioButton("Yes");
		rbtnOutGrpYes.setSelected(false);
		rbtnOutGrpYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutRootToggle(false);
			}
		});
		rbtnOutGrpYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnOutGrpYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnOutGrpYes, "cell 2 8");
		
		rbtnOutGrpNo = new JRadioButton("No");
		rbtnOutGrpNo.setSelected(true);
		rbtnOutGrpNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutRootToggle(true);
			}
		});
		rbtnOutGrpNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnOutGrpNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnOutGrpNo, "cell 2 8");
		
		lblOutRootNum = new JLabel("Type number of the outgroup:");
		lblOutRootNum.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutRootNum.setEnabled(false);
		panel.add(lblOutRootNum, "flowx,cell 2 9,alignx left");
		
		txtOutRootNum = new JTextField();
		txtOutRootNum.setColumns(5);
		txtOutRootNum.setText("1");
		txtOutRootNum.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutRootNum.setEnabled(false);
		panel.add(txtOutRootNum, "cell 2 9");
		
		lblGlobal = new JLabel("Global rearrangements:");
		lblGlobal.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGlobal, "flowx,cell 0 10 2 1,alignx right");
		
		rbtnGlobalYes = new JRadioButton("Yes");
		rbtnGlobalYes.setSelected(false);
		rbtnGlobalYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnGlobalYes.setBackground(new Color(204, 255, 255));
		rbtnGlobalYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalRearrToggle(true);
			}
		});
		panel.add(rbtnGlobalYes, "cell 2 10");
		
		rbtnGlobalNo = new JRadioButton("No");
		rbtnGlobalNo.setSelected(true);
		rbtnGlobalNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnGlobalNo.setBackground(new Color(204, 255, 255));
		rbtnGlobalNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GlobalRearrToggle(false);
			}
		});
		panel.add(rbtnGlobalNo, "cell 2 10");
		
		lblRandOrder = new JLabel("Randomize input order of species:");
		lblRandOrder.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOrder, "flowx,cell 0 11 2 1,alignx right");
		
		rbtnRandOrderYes = new JRadioButton("Yes");
		rbtnRandOrderYes.setSelected(false);
		rbtnRandOrderYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(true);
			}
		});
		rbtnRandOrderYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRandOrderYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnRandOrderYes, "cell 2 11");
		
		rbtnRandOrderNo = new JRadioButton("No, use input order");
		rbtnRandOrderNo.setSelected(true);
		rbtnRandOrderNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandOrderToggle(false);
			}
		});
		rbtnRandOrderNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRandOrderNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnRandOrderNo, "cell 2 11");
		
		lblRandSeed = new JLabel("Random number seed:");
		lblRandSeed.setFont(new Font("Arial", Font.BOLD, 13));
		lblRandSeed.setEnabled(false);
		panel.add(lblRandSeed, "flowx,cell 0 12 2 1,alignx right");
		
		txtRandSeed = new JTextField();
		txtRandSeed.setText("1");
		txtRandSeed.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRandSeed.setEnabled(false);
		txtRandSeed.setColumns(6);
		panel.add(txtRandSeed, "cell 2 12");
		
		lblRandOdd = new JLabel("(must be odd)");
		lblRandOdd.setEnabled(false);
		lblRandOdd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOdd, "cell 2 12");
		
		lblNumberJumble = new JLabel("Number of times to jumble:");
		lblNumberJumble.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberJumble.setEnabled(false);
		panel.add(lblNumberJumble, "flowx,cell 0 13 2 1,alignx right");

		txtNumberJumble = new JTextField();
		txtNumberJumble.setColumns(6);
		txtNumberJumble.setText("1");
		txtNumberJumble.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberJumble.setEnabled(false);
		panel.add(txtNumberJumble, "cell 2 13");
		
		lblMultSets = new JLabel("Analyze multiple data sets:");
		lblMultSets.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMultSets, "flowx,cell 0 14 2 1,alignx right");
		
		rbtnMultSetsYes = new JRadioButton("Yes");
		rbtnMultSetsYes.setSelected(false);
		rbtnMultSetsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AnalyzeMultToggle(true);
			}
		});
		rbtnMultSetsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultSetsYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnMultSetsYes, "cell 2 14");
		
		rbtnMultSetsNo = new JRadioButton("No");
		rbtnMultSetsNo.setSelected(true);
		rbtnMultSetsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AnalyzeMultToggle(false);
			}
		});
		rbtnMultSetsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultSetsNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnMultSetsNo, "cell 2 14");
		
		lblNumSet = new JLabel("Number of data sets:");
		lblNumSet.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumSet.setEnabled(false);
		panel.add(lblNumSet, "flowx,cell 0 15 2 1,alignx right");
		
		txtNumSets = new JTextField();
		txtNumSets.setColumns(6);
		txtNumSets.setText("1");
		txtNumSets.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumSets.setEnabled(false);
		panel.add(txtNumSets, "cell 2 15");
		
		lblPrintInput = new JLabel("Print out the data at start of run:");
		lblPrintInput.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInput, "flowx,cell 0 16 2 1,alignx right");
		
		rbtnPrintInputYes = new JRadioButton("Yes");
		rbtnPrintInputYes.setSelected(false);
		rbtnPrintInputYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		rbtnPrintInputYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintInputYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnPrintInputYes, "cell 2 16");
		
		rbtnPrintInputNo = new JRadioButton("No");
		rbtnPrintInputNo.setSelected(true);
		rbtnPrintInputNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rbtnPrintInputNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintInputNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnPrintInputNo, "cell 2 16");
		
		lblPrintTree = new JLabel("Print out tree:");
		lblPrintTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintTree, "flowx,cell 0 17 2 1,alignx right");
		
		rbtnPrintTreeYes = new JRadioButton("Yes");
		rbtnPrintTreeYes.setSelected(true);
		rbtnPrintTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(true);
			}
		});
		rbtnPrintTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintTreeYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnPrintTreeYes, "cell 2 17");
		
		rbtnPrintTreeNo = new JRadioButton("No");
		rbtnPrintTreeNo.setSelected(false);
		rbtnPrintTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(false);
			}
		});
		rbtnPrintTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintTreeNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnPrintTreeNo, "cell 2 17");
		
		lblTreeFile = new JLabel("Write out trees onto tree file:");
		lblTreeFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTreeFile, "flowx,cell 0 18 2 1,alignx right");
		
		rbtnTreeFileYes = new JRadioButton("Yes");
		rbtnTreeFileYes.setSelected(true);
		rbtnTreeFileYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(true);
			}
		});
		rbtnTreeFileYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnTreeFileYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnTreeFileYes, "cell 2 18");
		
		rbtnTreeFileNo = new JRadioButton("No");
		rbtnTreeFileNo.setSelected(false);
		rbtnTreeFileNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(false);
			}
		});
		rbtnTreeFileNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnTreeFileNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnTreeFileNo, "cell 2 18");

		
		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInd, "flowx,cell 0 19 2 1,alignx right");
		
		rbtnPrintIndYes = new JRadioButton("Yes");
		rbtnPrintIndYes.setBackground(new Color(204, 255, 255));
		rbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rbtnPrintIndYes.setSelected(true);
		panel.add(rbtnPrintIndYes, "cell 2 19");
		
		rbtnPrintIndNo = new JRadioButton("No");
		rbtnPrintIndNo.setBackground(new Color(204, 255, 255));
		rbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rbtnPrintIndNo.setSelected(false);
		panel.add(rbtnPrintIndNo, "cell 2 19");
		
		btnExecute = new JButton("Execute");
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				inputvals = new ContMLData();
				inputvals.infile = txtInFile.getText();
				inputvals.intree = txtInTree.getText();
				inputvals.outfile = txtOutFile.getText();
				inputvals.outfileopt = "w";
				inputvals.outtree = txtOutTree.getText();
				inputvals.outtreeopt = "w";
				inputvals.BestTree = rbtnSearchBestYes.isSelected();
				inputvals.UseLengths = rdbtnUseLengthsYes.isSelected();
				inputvals.GeneFreq = rbtnGeneFreq.isSelected();
				inputvals.AllAlleles = rbtnAllelYes.isSelected();
				inputvals.OutRoot = rbtnOutGrpYes.isSelected();
				inputvals.OutNum = Integer.parseInt(txtOutRootNum.getText());
				inputvals.GlobalRearr = rbtnGlobalYes.isSelected();
				inputvals.RandInput = rbtnRandOrderYes.isSelected();
				inputvals.RandNum = Integer.parseInt(txtRandSeed.getText());
				inputvals.Njumble = Integer.parseInt(txtNumberJumble.getText());
				inputvals.MultData = rbtnMultSetsYes.isSelected();
				inputvals.NumSets = Integer.parseInt(txtNumSets.getText());
				inputvals.PrintData = rbtnPrintInputYes.isSelected();
				inputvals.PrintInd = rbtnPrintIndYes.isSelected();
				inputvals.PrintTree = rbtnPrintTreeYes.isSelected();
				inputvals.WriteTree = rbtnTreeFileYes.isSelected();
				
				btnExecute.setEnabled(false);	
				String title = "Contml Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread contMLThread = new Thread() {
						public void run() {
							runContMLThreads();
						}
			  	    };
			  	    contMLThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		panel.add(btnExecute, "cell 2 20");
		
		btnQuit = new JButton("Quit");
		btnQuit.setHorizontalAlignment(SwingConstants.RIGHT);
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				if(phylipCall)
				{
					frmContMLControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 20,alignx left");
	}
	
	public boolean checkInputVals(){
		TestFileNames test = new TestFileNames();
		
		if (!test.DuplicateFileNames(inputvals.infile, "Input", inputvals.outfile, "Output"))
		{			
			return false;		
		}

		if (!test.FileAvailable(inputvals.infile, "Input"))
		{
			return false;
		}

		if (!inputvals.BestTree)
		{
			if (!test.FileAvailable(inputvals.intree, "Intree"))
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
		

		if (inputvals.WriteTree)
		{
			opt = test.FileAlreadyExists(inputvals.outtree, "Outtree");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inputvals.outfileopt = opt;
			}
		}
		
		if ((inputvals.RandNum % 2) == 0)
		{
			String msg1 = "Random number seed must be odd.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;			
		}
		return true;
	}
	
	protected void runContMLThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("contml", ContML.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("ContML");
    		return;
    	}
		try 
		{
	  	    Thread contMLRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			ContML contml = (ContML) Native.loadLibrary("contml", ContML.class);
		  			contml.contml(		
		  	    			inputvals.infile,
		  	    			inputvals.intree,
		  	    			inputvals.outfile,
		  	    			inputvals.outfileopt,
		  	    			inputvals.outtree,
		  	    			inputvals.outtreeopt,
		  	    			inputvals.BestTree,
		  	    			inputvals.UseLengths,
		  	    			inputvals.GeneFreq,
		  	    			inputvals.AllAlleles,
		  	       			inputvals.OutRoot,
		  	    			inputvals.OutNum,
		  	    			inputvals.GlobalRearr,
		  	    			inputvals.RandInput,
		  	    			inputvals.RandNum,
		  	    			inputvals.Njumble,
		  	     			inputvals.MultData,
		  	    			inputvals.NumSets,
		  	    			inputvals.PrintData,
		  	    			inputvals.PrintInd,
		  	    			inputvals.PrintTree,
		  	    			inputvals.WriteTree);
				    		
		  	    };
	  	    };
	  	    contMLRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (contMLRunThread.isAlive());
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
