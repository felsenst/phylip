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

public class KitschUserInterface {
	 public interface Kitsch extends Library {
	      public void kitsch(
			String infile,
			String intree,
			String outfile,
			String outfileopt,
			String outtree,
			String outtreeopt,
			String Method,
			boolean BestTree,
			double Power,
			boolean NegLengths,
			boolean LowerTMat,
			boolean UpperTMat,
			boolean Subreps,
			boolean RandInput,
			int RandNum,
			int Njumble,
			boolean MultData,
			int NumSets,
			boolean PrintData,
			boolean PrintInd,
			boolean PrintTree,
			boolean WriteTree
	 		);
	 }

	 public class KitschData {
		String infile;
		String intree;
		String outfile;
		String outfileopt;
		String outtree;
		String outtreeopt;
		String Method;
		boolean BestTree;
		double Power;
		boolean NegLengths;
		boolean LowerTMat;
		boolean UpperTMat;
		boolean Subreps;
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

	private KitschData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;

	private JFrame frmKitschControls;
	private JButton btnInputFile;
	private JTextField txtInputFile;
	private JTextField txtInputTree;
	private JButton btnInputTree;
	private JButton btnOutputFile;
	private JTextField txtOutputFile;
	private JButton btnOutputTree;
	private JTextField txtOutputTree;
	private JLabel lblMethod;
	private JComboBox cmbxMethod;
	private JLabel lblSearchBest;
	private JRadioButton rbtnSearchBestYes;
	private JRadioButton rbtnSearchBestNo;
	private JLabel lblMultSets;
	private JRadioButton rbtnMultSetsYes;
	private JRadioButton rbtnMultSetsNo;
	private JLabel lblLowTriMat;
	private JRadioButton rbtnLowTriMatYes;
	private JRadioButton rbtnLowTriMatNo;
	private JLabel lblUpTriMat;
	private JRadioButton rbtnUpTriMatYes;
	private JRadioButton rbtnUpTriMatNo;
	private JLabel lblSubreps;
	private JRadioButton rbtnSubrepsYes;
	private JRadioButton rbtnSubrepsNo;
	private JLabel lblRandOrder;
	private JRadioButton rbtnRandOrderYes;
	private JRadioButton rbtnRandOrderNo;
	private JLabel lblRandSeed;
	private JTextField txtRandSeed;
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
	private JLabel lblNumSet;
	private JTextField txtNumSets;
	private JLabel lblPower;
	private JTextField txtPower;
	private JLabel lblNegLen;
	private JRadioButton rbtnNegLenYes;
	private JRadioButton rbtnNegLenNo;
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
					KitschUserInterface window = new KitschUserInterface(args);
					window.frmKitschControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmKitschControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}
	
	protected void MethodToggle(int selected){
		if(selected == 0) //Fitch-Margoliash
		{	
			NegLenToggle(false);
		}
		else // Minimum Evolution
		{
			NegLenToggle(true);
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
			txtOutputTree.setEnabled(true);
		}
		else{
			rbtnTreeFileYes.setSelected(false);
			rbtnTreeFileNo.setSelected(true);
			btnOutputTree.setEnabled(false);
			txtOutputTree.setEnabled(false);
		}
	}
	
	protected void LowTriMatToggle(boolean isTrue){
		if(isTrue){
			rbtnLowTriMatYes.setSelected(true);
			rbtnLowTriMatNo.setSelected(false);
			if (rbtnUpTriMatYes.isSelected())
			{
				rbtnUpTriMatYes.setSelected(false);
				rbtnUpTriMatNo.setSelected(true);
			}
		}
		else{
			rbtnLowTriMatYes.setSelected(false);
			rbtnLowTriMatNo.setSelected(true);
		}
	}
	
	protected void UpTriMatToggle(boolean isTrue){
		if(isTrue){
			rbtnUpTriMatYes.setSelected(true);
			rbtnUpTriMatNo.setSelected(false);
			if (rbtnLowTriMatYes.isSelected())
			{
				rbtnLowTriMatYes.setSelected(false);
				rbtnLowTriMatNo.setSelected(true);
			}
		}
		else{
			rbtnUpTriMatYes.setSelected(false);
			rbtnUpTriMatNo.setSelected(true);
		}
	}
	
	protected void SubrepsToggle(boolean isTrue){
		if(isTrue){
			rbtnSubrepsYes.setSelected(true);
			rbtnSubrepsNo.setSelected(false);
		}
		else{
			rbtnSubrepsYes.setSelected(false);
			rbtnSubrepsNo.setSelected(true);
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
	
	protected void IntreeToggle(boolean useIntree){
		if(useIntree){
			txtInputTree.setEnabled(true);
			btnInputTree.setEnabled(true);
			rbtnSearchBestYes.setSelected(false);
			rbtnSearchBestNo.setSelected(true);
		}
		else{
			txtInputTree.setEnabled(false);
			btnInputTree.setEnabled(false);
			rbtnSearchBestYes.setSelected(true);
			rbtnSearchBestNo.setSelected(false);
		}
	}
	
	protected void NegLenToggle(boolean isAllowed) {
		if (isAllowed) {
			rbtnNegLenYes.setSelected(true);
			rbtnNegLenNo.setSelected(false);
		} else {
			rbtnNegLenYes.setSelected(false);
			rbtnNegLenNo.setSelected(true);
		}
	}

	/**
	 * Create the application.
	 */
	public KitschUserInterface(String[] args) {
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

		frmKitschControls = new JFrame();
		frmKitschControls.setBackground(new Color(204, 255, 255));
		frmKitschControls.setTitle("Kitsch");
		frmKitschControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmKitschControls.setBounds(100, 100, 625, 680);
		frmKitschControls.setPreferredSize(new Dimension(frmKitschControls.getBounds().width, frmKitschControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmKitschControls.getPreferredSize());
		frmKitschControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmKitschControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow]", "[][][][]"));
		
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
		txtInputFile.setBounds(168, 11, 412, 20);
		panel.add(txtInputFile, "cell 1 0 2 1,growx");
		
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
		txtInputTree.setBounds(168, 36, 412, 20);
		panel.add(txtInputTree, "cell 1 1 2 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		panel.add(btnOutputFile, "cell 0 2,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 2 2 1,growx");
		

		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputTree);
			}
		});
		panel.add(btnOutputTree, "cell 0 3,growx");
		
		txtOutputTree = new JTextField();
		txtOutputTree.setText("outtree");
		txtOutputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutputTree, "cell 1 3 2 1,growx");
		
		lblMethod = new JLabel("Method:");
		lblMethod.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMethod.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMethod, "flowx,cell 0 4 2 1,alignx right");
		
		cmbxMethod = new JComboBox();
		cmbxMethod.setModel(new DefaultComboBoxModel(new String[] {"Fitch-Margoliash", "Minimum Evolution"}));
		cmbxMethod.setSelectedIndex(0);
		cmbxMethod.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				MethodToggle(cmbxMethod.getSelectedIndex());
			}
		});
		cmbxMethod.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxMethod, "cell 2 4,growx");
		
		lblSearchBest = new JLabel("Search for best tree:");
		lblSearchBest.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSearchBest.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSearchBest, "flowx,cell 0 5 2 1,alignx right");

		rbtnSearchBestYes = new JRadioButton("Yes");
		rbtnSearchBestYes.setSelected(true);
		rbtnSearchBestYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSearchBestYes.setBackground(new Color(204, 255, 255));
		rbtnSearchBestYes.setBounds(274, 146, 56, 23);
		rbtnSearchBestYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(false);
			}
		});
		panel.add(rbtnSearchBestYes, "cell 2 5");
		
		rbtnSearchBestNo = new JRadioButton("No, use user input tree");
		rbtnSearchBestNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSearchBestNo.setBackground(new Color(204, 255, 255));
		rbtnSearchBestNo.setBounds(329, 146, 220, 23);
		rbtnSearchBestNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(true);
			}
		});
		panel.add(rbtnSearchBestNo, "cell 2 5");
		
		lblPower = new JLabel("Power:");
		lblPower.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPower.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPower, "flowx,cell 0 6 2 1,alignx right");
		
		txtPower = new JTextField();
		txtPower.setText("2.0");
		txtPower.setFont(new Font("Arial", Font.PLAIN, 13));
		txtPower.setColumns(6);
		panel.add(txtPower, "cell 2 6");
		
		lblNegLen = new JLabel("Negative branch lengths allowed:");
		lblNegLen.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNegLen.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNegLen, "flowx,cell 0 7 2 1,alignx right");
	
		rbtnNegLenYes = new JRadioButton("Yes");
		rbtnNegLenYes.setSelected(false);
		rbtnNegLenYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnNegLenYes.setBackground(new Color(204, 255, 255));
		rbtnNegLenYes.setBounds(274, 196, 56, 23);
		rbtnNegLenYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				NegLenToggle(true);
			}
		});
		panel.add(rbtnNegLenYes, "cell 2 7");
		
		rbtnNegLenNo = new JRadioButton("No");
		rbtnNegLenNo.setSelected(true);
		rbtnNegLenNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnNegLenNo.setBackground(new Color(204, 255, 255));
		rbtnNegLenNo.setBounds(329, 196, 50, 23);
		rbtnNegLenNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				NegLenToggle(false);
			}
		});
		panel.add(rbtnNegLenNo, "cell 2 7");
		
		lblLowTriMat = new JLabel("Lower-triangular data matrix:");
		lblLowTriMat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblLowTriMat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblLowTriMat, "flowx,cell 0 8 2 1,alignx right");
		
		rbtnLowTriMatYes = new JRadioButton("Yes");
		rbtnLowTriMatYes.setSelected(false);
		rbtnLowTriMatYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				LowTriMatToggle(true);
			}
		});
		rbtnLowTriMatYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnLowTriMatYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnLowTriMatYes, "cell 2 8");
		
		rbtnLowTriMatNo = new JRadioButton("No");
		rbtnLowTriMatNo.setSelected(true);
		rbtnLowTriMatNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				LowTriMatToggle(false);
			}
		});
		rbtnLowTriMatNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnLowTriMatNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnLowTriMatNo, "cell 2 8");
		
		lblUpTriMat = new JLabel("Upper-triangular data matrix:");
		lblUpTriMat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblUpTriMat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUpTriMat, "flowx,cell 0 9 2 1,alignx right");
		
		rbtnUpTriMatYes = new JRadioButton("Yes");
		rbtnUpTriMatYes.setSelected(false);
		rbtnUpTriMatYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UpTriMatToggle(true);
			}
		});
		rbtnUpTriMatYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnUpTriMatYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnUpTriMatYes, "cell 2 9");
		
		rbtnUpTriMatNo = new JRadioButton("No");
		rbtnUpTriMatNo.setSelected(true);
		rbtnUpTriMatNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				UpTriMatToggle(false);
			}
		});
		rbtnUpTriMatNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnUpTriMatNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnUpTriMatNo, "cell 2 9");
		
		lblSubreps = new JLabel("Subreplicates:");
		lblSubreps.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSubreps.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSubreps, "flowx,cell 0 10 2 1,alignx right");
		
		rbtnSubrepsYes = new JRadioButton("Yes");
		rbtnSubrepsYes.setSelected(false);
		rbtnSubrepsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SubrepsToggle(true);
			}
		});
		rbtnSubrepsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSubrepsYes.setBackground(new Color(204, 255, 255));
		panel.add(rbtnSubrepsYes, "cell 2 10");
		
		rbtnSubrepsNo = new JRadioButton("No");
		rbtnSubrepsNo.setSelected(true);
		rbtnSubrepsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SubrepsToggle(false);
			}
		});
		rbtnSubrepsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSubrepsNo.setBackground(new Color(204, 255, 255));
		panel.add(rbtnSubrepsNo, "cell 2 10");
		
		lblRandOrder = new JLabel("Randomize input order of species:");
		lblRandOrder.setHorizontalAlignment(SwingConstants.RIGHT);
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
		lblRandSeed.setHorizontalAlignment(SwingConstants.RIGHT);
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
		lblRandOdd.setHorizontalAlignment(SwingConstants.LEFT);
		lblRandOdd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOdd, "cell 2 12");
		
		lblNumberJumble = new JLabel("Number of times to jumble:");
		lblNumberJumble.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberJumble.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumberJumble.setEnabled(false);
		panel.add(lblNumberJumble, "flowx,cell 0 13 2 1,alignx right");

		txtNumberJumble = new JTextField();
		txtNumberJumble.setColumns(6);
		txtNumberJumble.setText("1");
		txtNumberJumble.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberJumble.setEnabled(false);
		txtNumberJumble.setBounds(274, 351, 86, 20);
		panel.add(txtNumberJumble, "cell 2 13");

		lblMultSets = new JLabel("Analyze multiple data sets:");
		lblMultSets.setHorizontalAlignment(SwingConstants.RIGHT);
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
		lblNumSet.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumSet.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumSet.setEnabled(false);
		panel.add(lblNumSet, "flowx,cell 0 15 2 1,alignx right");
		
		txtNumSets = new JTextField();
		txtNumSets.setText("1");
		txtNumSets.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumSets.setEnabled(false);
		txtNumSets.setColumns(6);
		panel.add(txtNumSets, "cell 2 15");

		lblPrintInput = new JLabel("Print out the data at start of run:");
		lblPrintInput.setHorizontalAlignment(SwingConstants.RIGHT);
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
		lblPrintTree.setHorizontalAlignment(SwingConstants.RIGHT);
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
		lblTreeFile.setHorizontalAlignment(SwingConstants.RIGHT);
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
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
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

				inputvals = new KitschData();
				inputvals.infile = txtInputFile.getText();
				inputvals.intree = txtInputTree.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.outtree = txtOutputTree.getText();
				inputvals.outtreeopt = "w";
				if (cmbxMethod.getSelectedIndex() == 0)
				{
					inputvals.Method = "FM";
				}
				else
				{
					inputvals.Method = "ME";				
				}
				inputvals.BestTree = rbtnSearchBestYes.isSelected();
				inputvals.Power = Double.parseDouble(txtPower.getText());
				inputvals.NegLengths = rbtnNegLenYes.isSelected();
				inputvals.LowerTMat = rbtnLowTriMatYes.isSelected();
				inputvals.UpperTMat = rbtnUpTriMatYes.isSelected();
				inputvals.Subreps   = rbtnSubrepsYes.isSelected();
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
				String title = "Kitsch Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread kitschThread = new Thread() {
						public void run() {
							runKitschThreads();
						}
			  	    };
			  	  kitschThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		panel.add(btnExecute, "cell 2 20,alignx center");
		
		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				if(phylipCall)
				{
					frmKitschControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 20");
	}
	
	public boolean checkInputVals(){
		
		// check files
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
		
		// check data
		if (inputvals.LowerTMat && inputvals.UpperTMat)
		{
			String msg1 = "Data matrix cannot be both upper and lower triangular.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if ((inputvals.RandNum % 2) == 0)
		{
			String msg1 = "Random number seed must be odd.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;			
		}
		return true;
		
	 }
	
	protected void runKitschThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("kitsch", Kitsch.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("Kitsch");
    		return;
    	}
		try 
		{
	  	    Thread kitschRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			Kitsch kitsch = (Kitsch) Native.loadLibrary("kitsch", Kitsch.class);
		  	        kitsch.kitsch(		
		  	    			inputvals.infile,
		  	    			inputvals.intree,
		  	    			inputvals.outfile,
		  	    			inputvals.outfileopt,
		  	    			inputvals.outtree,
		  	    			inputvals.outtreeopt,
		  	    			inputvals.Method,
		  	    			inputvals.BestTree,
		  	     			inputvals.Power,
		  	    			inputvals.NegLengths,
		  	    			inputvals.LowerTMat,
		  	    			inputvals.UpperTMat,
		  	    			inputvals.Subreps,
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
	  	    kitschRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (kitschRunThread.isAlive());
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
