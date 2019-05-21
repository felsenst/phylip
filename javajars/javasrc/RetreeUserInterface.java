package phylip;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.print.Book;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterJob;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Scanner;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ActionMap;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

import utilities.TestFileNames;

import dnaml.DnaMLUserInterface;
import drawgram.DrawgramUserInterface;
import drawtree.DrawtreeUserInterface;

@SuppressWarnings("serial")
public class RetreeUserInterface {
	
	private JFrame frmRetreeControls;
	private JPanel displayPanel;
	private JMenuBar menuFile;
	
	private JMenu mnTree;
	private JMenuItem mntmRead;
	private JMenuItem mntmSave;
	private JMenuItem mntmPrint;
	private JMenuItem mntmAbout;
	private JMenuItem mntmQuit;
	private JMenuItem mntmUndo;
	private JMenuItem mntmClade;
	
	private JMenu mnPhylip;
	private JMenuItem mntmDnaML;
	private JMenuItem mntmDrawGram;
	private JMenuItem mntmDrawTree;
	
	private JPopupMenu popAction;
	private JMenuItem pntmMoveNode;
	private JMenuItem pntmOutgroup;
	private JMenuItem pntmMidpoint;
	private JMenuItem pntmTranspose;
	private JMenuItem pntmFlip;
	private JMenuItem pntmLength;
	private JMenuItem pntmShowLen;
	private JMenuItem pntmName;
	private JMenuItem pntmUndo;
	
	private JPopupMenu popMove;
	private JMenuItem mntmMoveTo;
	private JMenuItem mntmMoveBefore;
	private JMenuItem mntmNoMove;
	
	private JFrame frmRename;
	private JLabel lblRename;
	private JTextField txtRename;
	private JButton btnRenameOK;
	private JButton btnRenameCancel;
	
	private JFrame frmNewLen;
	private JLabel lblNewLen;
	private JTextField txtNewLen;
	private JButton btnNewLenOK;
	private JButton btnNewLenCancel;

	private JFrame frmAbout;
	private JPanel contentPane;
	private JButton btnOK;
	
	private String filedir;
	private boolean phylipCall;
	
	// Newick tree info
	int m_curr_char;
	char[] m_newick;
	NewickNode root;
	NewickNode newNode;
	NewickNode newBase;
	public enum NodeType{ROOT, NODE, TIP, UNDEF}
	int nodeNum;
	boolean lengthFound;
    String outputTree;
    String outputTreePath;

	// plot parameters
	Boolean treeRead;
	int numTips; 
	int maxNameLen;
	ArrayList<NewickNode> treeTips;
	int xOffset;
	int rootLen;
    double xScaleFactor;
    double yScaleFactor;
    int nameOffsetX;
    int nameTweakX;
    int nameTweakY;
    ArrayList<BranchSeg> branchSegs;
    BranchSeg newseg;
    double tipYcor;
    Boolean showLens;
    
    // mouse parameters
    Point2D.Float mousePt;
	int bestTip;
	int bestSeg;
	int moveSeg;
	int moveEnd;
	Component mousePick;
	
	// edit control
	Boolean doMove;
	Cursor oldCursor;
    NewickNode editNode;
	NewickNode longestBranch;
	NewickNode nextLongBranch;

	// Undo array
	ArrayList<String> UndoTrees;
	
	// keyboard shortcuts
	KeyStroke ctlA;
	KeyStroke ctlP;
	KeyStroke ctlQ;
	KeyStroke ctlR;
	KeyStroke ctlS;
	KeyStroke ctlZ;

	Action printTree;
	Action readTree;
	Action saveTree;
	Action quitRun;
	Action showAbout;
	Action doUndo;

	
	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					RetreeUserInterface window = new RetreeUserInterface(args);
					window.frmRetreeControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	/**
	 * Create the application.
	 */
	public RetreeUserInterface(String[] args) {
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
		// initialize data
		final String[] phylip = new String[] {"Phylip"};

		filedir = System.getProperty("user.dir");
		treeRead = false;
		doMove = false;
		showLens = false;
	    branchSegs = new ArrayList<BranchSeg>();
		UndoTrees = new ArrayList<String> ();	
   		outputTree = "RetreeOuttree";
   		outputTreePath = "RetreeOuttree.txt";

		frmRetreeControls = new JFrame();
        frmRetreeControls.setTitle("Retree");
        frmRetreeControls.setBackground(new Color(204, 255, 255));
        frmRetreeControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frmRetreeControls.getContentPane().setLayout(null);
        frmRetreeControls.setVisible(true);

        frmRetreeControls.setBounds(100, 100, 600, 600);
   	 	contentPane = new JPanel();
   	 	contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
   	 	contentPane.setLayout(null);
   	 	frmRetreeControls.add(contentPane);
   	 	
   	 	displayPanel = new JPanel()
        {
            public void paintComponent(Graphics graph)
            {
                draw(graph);
            }
        };
	    displayPanel.setBackground(Color.WHITE);
	    displayPanel.setBounds(0, 0, 5000, 5000);
	    displayPanel.setLayout(null);
   	 	frmRetreeControls.add(displayPanel);
	    
		// define accelerator keys
	    // ***************************************************************************
	    // Note: the accelerator keys are registered twice:
	    // Once with the Action Map (so they work anywhere on the screen at any time)
	    // Once with the MenuItem they pertain to
	    // It works fine with only the Action Map, but it doesn't show the accelerator 
	    // key combination in the MenuItem unless setAccelerator is called.  
	    // JRM
	    // ***************************************************************************

	    int shortcut = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(); // make key definitions platform independent
		
		ctlA = KeyStroke.getKeyStroke(KeyEvent.VK_A,shortcut); // About Retree
		ctlP = KeyStroke.getKeyStroke(KeyEvent.VK_P,shortcut); // Print tree
		ctlQ = KeyStroke.getKeyStroke(KeyEvent.VK_Q,shortcut); // Quit
		ctlR = KeyStroke.getKeyStroke(KeyEvent.VK_R,shortcut); // Read tree
		ctlS = KeyStroke.getKeyStroke(KeyEvent.VK_S,shortcut); // Save tree
		ctlZ = KeyStroke.getKeyStroke(KeyEvent.VK_Z,shortcut); // Undo
	    
		showAbout = new ShowAbout();
		printTree = new PrintTree();
		quitRun   = new QuitRun();
		readTree  = new ReadTree();
		saveTree  = new SaveTree();
		doUndo    = new UndoLast();
		
		registerShortcuts();
		
		oldCursor = displayPanel.getCursor();
	    displayPanel.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e){
				if(doMove)
				{
					popMove.show(e.getComponent(), e.getX(), e.getY());
					mousePt = new Point2D.Float(e.getX(), e.getY());
					mousePick = e.getComponent();
					displayPanel.setCursor(oldCursor);
					double minsegdist = 100000;
					moveSeg = 0;
					for(int i=0; i<branchSegs.size(); i++)
					{
						BranchSeg seg = branchSegs.get(i);
						if ((seg.distanceToBranch(mousePt)) < minsegdist)
						{
							minsegdist = seg.distanceToBranch(mousePt);
							moveSeg = i;
						}
					}		
			        branchSegs.get(moveSeg).m_endnode.m_color = Color.GREEN;
					doMove = false;
					frmRetreeControls.repaint();        					
				}
				else
				{
					if (root != null)
						{
		 				root.ClearTree(root);
						popAction.show(e.getComponent(), e.getX(), e.getY());
						mousePt = new Point2D.Float(e.getX(), e.getY());
						//System.out.println("Mouse Clicked at X: " + e.getX() + " Y: " + e.getY());
						double mintipdist = 100000; 
						bestTip = 0;
						for(int i=0; i<treeTips.size(); i++)
						{
							NewickNode node = treeTips.get(i);
							if (node.m_dispPt.distance(mousePt) < mintipdist)
							{
								mintipdist = node.m_dispPt.distance(mousePt);
								bestTip = i;
							}
						}
						
						double minsegdist = 100000;
						bestSeg = 0;
						for(int i=0; i<branchSegs.size(); i++)
						{
							BranchSeg seg = branchSegs.get(i);
							if ((seg.distanceToBranch(mousePt)) < minsegdist)
							{
								minsegdist = seg.distanceToBranch(mousePt);
								bestSeg = i;
							}
						}
						
						if (mintipdist <= minsegdist)
						{
					        pntmName.setEnabled(true);
					        pntmLength.setEnabled(false);
					        pntmShowLen.setEnabled(true);
							pntmFlip.setEnabled(false);
							pntmTranspose.setEnabled(false);
					        pntmMoveNode.setEnabled(false);
					        pntmOutgroup.setEnabled(false);
					        mntmClade.setEnabled(false);
					        pntmMidpoint.setEnabled(true);
					        treeTips.get(bestTip).m_color = Color.BLUE;
						}
						else
						{
					        pntmName.setEnabled(false);
					        pntmLength.setEnabled(true);
					        pntmShowLen.setEnabled(true);
					        pntmOutgroup.setEnabled(true);
					        mntmClade.setEnabled(false);
					        pntmMidpoint.setEnabled(true);
					        
					        branchSegs.get(bestSeg).m_endnode.m_color = Color.BLUE;
					        if (branchSegs.get(bestSeg).m_endnode.m_type != NodeType.TIP)
					        {
								pntmFlip.setEnabled(true);
								pntmTranspose.setEnabled(true);
					        }
					        else
					        {
								pntmFlip.setEnabled(false);
								pntmTranspose.setEnabled(false);
					        }
					        
					        if (branchSegs.get(bestSeg).m_endnode.m_type == NodeType.ROOT)
					        {
						        pntmMoveNode.setEnabled(false);				        	
					        }
					        else
					        {
						        pntmMoveNode.setEnabled(true);
					        }
						}
						frmRetreeControls.repaint();        					
					}	
				}
			}
	    });

	    frmRetreeControls.getContentPane().add(displayPanel);
 	 	         
        menuFile = new JMenuBar();
        frmRetreeControls.setJMenuBar(menuFile);
        
        mnTree = new JMenu("Tree");
        menuFile.add(mnTree);
         
        mntmRead = new JMenuItem("Read");
        mntmRead.setEnabled(true);
        mntmRead.setAccelerator(ctlR);
        mntmRead.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		UndoTrees.clear();
 				readTree();
		        mntmSave.setEnabled(true);
		        mntmPrint.setEnabled(true);
		        pntmName.setEnabled(true);
		        mntmDnaML.setEnabled(true);
		        mntmDrawGram.setEnabled(true);
		        mntmDrawTree.setEnabled(true);
				frmRetreeControls.repaint();        					
		        doMove = false;
      	}
       });
        mnTree.add(mntmRead);
       
        mntmSave = new JMenuItem("Save");
        mntmSave.setEnabled(false);
        mntmSave.setAccelerator(ctlS);
        mntmSave.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		saveTree();
				doMove = false;
       	}
       });
        mnTree.add(mntmSave);

        mntmPrint = new JMenuItem("Print tree");
        mntmPrint.setEnabled(false);
        mntmPrint.setAccelerator(ctlP);
        mntmPrint.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
		        printTree();
				doMove = false;
		    }
		});      
        mnTree.add(mntmPrint);        

        mntmUndo = new JMenuItem("Undo");
        mntmUndo.setEnabled(false);
        mntmUndo.setAccelerator(ctlZ);
        mntmUndo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				undoLast();
		    }
		});
        mnTree.add(mntmUndo);     

        mntmAbout = new JMenuItem("About");
        mntmAbout.setAccelerator(ctlA);
        mntmAbout.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		showAbout();
        	}
        });
        mnTree.add(mntmAbout);
        
        frmRetreeControls.setVisible(true);

        mntmQuit = new JMenuItem("Quit");
        mntmQuit.setAccelerator(ctlQ);
        mntmQuit.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmRetreeControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
            }
        });
        
        mnTree.add(mntmQuit);
        
        mnPhylip = new JMenu("Phylip");
        menuFile.add(mnPhylip);
        
        mntmDnaML = new JMenuItem("Dnaml");
        mntmDnaML.setEnabled(false);
        mntmDnaML.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		saveInterumTree();		        
				try {
					DnaMLUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
			}
        });
        mnPhylip.add(mntmDnaML);
        
        mntmDrawGram = new JMenuItem("Drawgram");
        mntmDrawGram.setEnabled(false);
        mntmDrawGram.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		saveInterumTree();
				try {
					DrawgramUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
			}
        });
        mnPhylip.add(mntmDrawGram);
       
        mntmDrawTree = new JMenuItem("Drawtree");
        mntmDrawTree.setEnabled(false);
        mntmDrawTree.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		saveInterumTree();
				try {
					DrawtreeUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
    		}
       });
        mnPhylip.add(mntmDrawTree);
    
        popAction= new JPopupMenu();
        
        pntmMoveNode = new JMenuItem("Move node");
        pntmMoveNode.setEnabled(false);
        pntmMoveNode.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e){
 				root.ClearTree(root);
				frmRetreeControls.repaint();
			}
		});
        pntmMoveNode.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				doMove = true;
				displayPanel.setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
		    }
		});
        
        popAction.add(pntmMoveNode);

        popMove= new JPopupMenu();

        mntmMoveBefore = new JMenuItem("Before node");
        mntmMoveBefore.setEnabled(true);
        mntmMoveBefore.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
	       		NewickNode movingNode = branchSegs.get(bestSeg).m_endnode;
	       		NewickNode moveBeforeNode = branchSegs.get(moveSeg).m_endnode;
	       		moveNodeBefore(movingNode, moveBeforeNode);
				Redisplay(true);
		    }
		});        
        popMove.add(mntmMoveBefore);

        mntmMoveTo = new JMenuItem("To node");
        mntmMoveTo.setEnabled(true);
        mntmMoveTo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
	       		NewickNode movingNode = branchSegs.get(bestSeg).m_endnode;
	       		NewickNode moveToNode;
	       		if ((Math.abs(branchSegs.get(moveSeg).m_startnode.m_dispPt.x-mousePt.x)) < 
	       			 Math.abs((branchSegs.get(moveSeg).m_endnode.m_dispPt.x-mousePt.x)))
       			{
	       			moveToNode = branchSegs.get(moveSeg).m_startnode;
       			}
	       		else if(branchSegs.get(moveSeg).m_endnode.m_type == NodeType.TIP)
	       		{
	       			moveToNode = branchSegs.get(moveSeg).m_startnode;	       			
	       		}
	       		else
	       		{
	       			moveToNode = branchSegs.get(moveSeg).m_endnode;	       				       			
	       		}
	       		moveNodeTo(movingNode, moveToNode);
  				calcNodeCoords(root);
				Redisplay(true);
		    }
		});        
        popMove.add(mntmMoveTo);

        mntmNoMove = new JMenuItem("Cancel");
        mntmNoMove.setEnabled(true);
        mntmNoMove.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				root.ClearTree(root);
				frmRetreeControls.repaint();
		    }
		});        
        popMove.add(mntmNoMove);

        pntmOutgroup = new JMenuItem("Outgroup");
        pntmOutgroup.setEnabled(false);
        pntmOutgroup.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				doMove = false;
	       		NewickNode movingNode = branchSegs.get(bestSeg).m_endnode;
	       		outgroupRoot(movingNode);
				Redisplay(true);
		    }
		});
        popAction.add(pntmOutgroup);

        pntmMidpoint = new JMenuItem("Midpoint root");
        pntmMidpoint.setEnabled(false);
        pntmMidpoint.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				midpointRoot();
				Redisplay(true);
		    }
		});
        popAction.add(pntmMidpoint);

        pntmTranspose = new JMenuItem("Transpose");
        pntmTranspose.setEnabled(false);
        pntmTranspose.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
	       		editNode = branchSegs.get(bestSeg).m_endnode;
	       		transposeBranch(editNode);
				Redisplay(true);
		    }
		});
        popAction.add(pntmTranspose);

        pntmFlip = new JMenuItem("Flip");
        pntmFlip.setEnabled(false);
        pntmFlip.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
	       		editNode = branchSegs.get(bestSeg).m_endnode;
	       		flipBranch(editNode);
				Redisplay(true);
		    }
		});
        
        popAction.add(pntmFlip);

        pntmLength = new JMenuItem("Edit branch length");
        pntmLength.setEnabled(false);
        pntmLength.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
        		frmNewLen = new JFrame();
        		frmNewLen.setVisible(true);
        		frmNewLen.setBounds(20, 50, 300, 130);
        		frmNewLen.setTitle("Edit Name");
        		frmNewLen.setLayout(null);
        		frmNewLen.setResizable(false);
        		
        		lblNewLen = new JLabel("Enter new length");
        		lblNewLen.setBounds(50, 10, 200, 14);
        		frmNewLen.add(lblNewLen);
        		
        		txtNewLen = new JTextField();
        		editNode = branchSegs.get(bestSeg).m_endnode;
        		txtNewLen.setText(editNode.m_length.toString());
        		txtNewLen.setBounds(35, 45, 200, 20);
        		frmNewLen.add(txtNewLen);
        		
        		btnNewLenOK = new JButton();
        		btnNewLenOK.setText("OK");
        		btnNewLenOK.setBounds(205, 75, 84, 25);
        		btnNewLenOK.addActionListener(new ActionListener() {
        			public void actionPerformed(ActionEvent e) {
        				lengthFound = true;
        	       		editNode.m_length =  Double.parseDouble(txtNewLen.getText());
          				frmNewLen.dispose();
          				calcNodeCoords(root);
          				//root.DumpTree(root);
        				Redisplay(true);
        			}
        		});
        		frmNewLen.add(btnNewLenOK);
        		
        		btnNewLenCancel = new JButton();
        		btnNewLenCancel.setText("Cancel");
        		btnNewLenCancel.setBounds(125, 75, 84, 25);
        		btnNewLenCancel.addActionListener(new ActionListener() {
        			public void actionPerformed(ActionEvent e) {
           				frmNewLen.dispose();
          				root.ClearTree(root);
     					frmRetreeControls.repaint();
      			    }
        		});
        		frmNewLen.add(btnNewLenCancel);	
			}
		});
        
        popAction.add(pntmLength);
        
        pntmShowLen = new JMenuItem("Show lengths");
        pntmShowLen.setEnabled(false);
        pntmShowLen.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (showLens)
				{
					showLens = false;
				}
				else
				{
					showLens = true;
				}
   				Redisplay(false);
		    }
		});
	        
	    popAction.add(pntmShowLen);

        pntmName = new JMenuItem("Edit name");
        pntmName.setEnabled(false);
        pntmName.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				//System.out.println("Edit Name clicked for item near x: " + mouseX + " y: " + mouseY);
				
				frmRename = new JFrame();
        		frmRename.setVisible(true);
        		frmRename.setBounds(20, 50, 300, 130);
        		frmRename.setTitle("Edit Name");
        		frmRename.setLayout(null);
        		frmRename.setResizable(false);
        		
        		lblRename = new JLabel("Enter new name");
        		lblRename.setBounds(50, 10, 200, 14);
        		frmRename.add(lblRename);
        		
        		txtRename = new JTextField();
        		editNode = treeTips.get(bestTip);
        		txtRename.setText(editNode.m_name);
        		txtRename.setBounds(35, 45, 200, 20);
        		frmRename.add(txtRename);
        		
        		btnRenameOK = new JButton();
        		btnRenameOK.setText("OK");
        		btnRenameOK.setBounds(205, 75, 84, 25);
        		btnRenameOK.addActionListener(new ActionListener() {
        			public void actionPerformed(ActionEvent e) {
        	       		editNode.m_name = txtRename.getText();
          				frmRename.dispose();
        				Redisplay(true);
        			}
        		});
        		frmRename.add(btnRenameOK);
        		
        		btnRenameCancel = new JButton();
        		btnRenameCancel.setText("Cancel");
        		btnRenameCancel.setBounds(125, 75, 84, 25);
        		btnRenameCancel.addActionListener(new ActionListener() {
        			public void actionPerformed(ActionEvent e) {
           				frmRename.dispose();
          				root.ClearTree(root);
     					frmRetreeControls.repaint();
      			    }
        		});
        		frmRename.add(btnRenameCancel);
 			
 				
		    }
		});
        popAction.add(pntmName);
        
        mntmClade = new JMenuItem("Clade");
        mntmClade.setEnabled(false);
        mntmClade.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				doMove = false;
		    }
		});
        //popAction.add(mntmClade);

        pntmUndo = new JMenuItem("Undo");
        pntmUndo.setEnabled(false);
        pntmUndo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				undoLast();
		    }
		});
        popAction.add(pntmUndo);     
	}
	
	void draw(Graphics g) {
	    Graphics2D g2d = (Graphics2D) g;
		g2d.setColor(frmRetreeControls.getBackground());
		g2d.fillRect(0, 0, frmRetreeControls.getWidth(), frmRetreeControls.getHeight());
		g2d.setColor(Color.BLACK);

		plotRetree(g2d);
	}
	
	void plotRetree(Graphics2D g2d){
		//BasicStroke stroke1 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] {4.0f, 4.0f}, 0.0f);
		//BasicStroke stroke2 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 1.0f);
		BasicStroke stroke3 = new BasicStroke(5.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 3.0f);
	    Font plotFont = new Font("helvetica", Font.BOLD, 14);
	    g2d.setFont(plotFont);
		g2d.setStroke(stroke3);
		
		//System.out.println("In plotRetree");
	    
	    if (treeRead)
	    {
	    	xOffset = 20;
	    	rootLen = 30;
	    	nameOffsetX = 5;
	    	nameTweakX = 30;
	    	nameTweakY = 5;
		    			
			// find longest path and longest label
			FontMetrics fm = g2d.getFontMetrics(plotFont);
			int nameWid = 0;
			double pathWid = 0.0;
			for(int i=0; i<treeTips.size(); i++)
			{
				if (fm.stringWidth(treeTips.get(i).m_name) > nameWid)
				{
					nameWid = fm.stringWidth(treeTips.get(i).m_name);
				}
				if (treeTips.get(i).m_corPt.x > pathWid)
				{
					pathWid = treeTips.get(i).m_corPt.x;
				}
			}
			
		    xScaleFactor = (frmRetreeControls.getWidth() - (nameWid + xOffset + rootLen + nameTweakX)) / pathWid;
		    yScaleFactor = frmRetreeControls.getHeight() / (double)(treeTips.size() + 1);
						
			// plot tree
		    if (!branchSegs.isEmpty())
		    {
		    	branchSegs.clear();
		    }
			root.PlotTree(root, g2d);
	    }
	}
   
	public void readTree() {
	    JFileChooser fileChooser = new JFileChooser(filedir);

		int option = fileChooser.showOpenDialog(frmRetreeControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			String intree = selectedFile.getPath();
			String newickTree = "";
			nodeNum = -1;
		    try 
		    {
				Scanner sf = new Scanner(new File(intree));
				while (sf.hasNextLine()) 
				{
					String curline = sf.nextLine();
					newickTree += curline;
				}
				addToUndo(newickTree);
				parseNewick(newickTree);
				String title = "Retree: " + selectedFile.getName();
		        frmRetreeControls.setTitle(title);
		    }
			catch (FileNotFoundException e)
			{
				String msg = "Input tree: ";
				msg += intree;
				msg += " does not exist.";
				JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);			
			}
		}	
	}
	
	public void parseNewick(String tree){
		m_newick = new char[tree.length()];
		tree.getChars(0, tree.length(), m_newick, 0);
		m_curr_char = 0;
		root = new NewickNode();
		root.m_type = NodeType.ROOT;
	    nodeNum = 0;
	    root.m_nodenum = nodeNum;
		NewickNode current = root;
	    treeTips = new ArrayList<NewickNode>();	
	    lengthFound = false; // if even one is found on input, the whole output tree gets lengths
		
		// parameters needed for plot
		numTips = 0;
		maxNameLen = 0;
		treeRead = false;
	
		while (m_newick[m_curr_char] != ';')
		{
			switch(m_newick[m_curr_char])
			{
            case '(':
                ++m_curr_char;
                // create first daughter
                current = current.AddChild();
                //System.out.print("found first daughter: ");
                break;
            case ',':
                ++m_curr_char;
                // create additional daughter
                current = current.GetParent().AddChild();
                //System.out.print("found second daughter: ");
                break;
            case ')':
                ++m_curr_char;
                // coalesce daughters
                current = current.GetParent().Terminate();
                //System.out.println("end of pair");
                break;
            case ' ':
            case '\n':
            case '\t':
            case '\r':
                ++m_curr_char;
                // skip whitespace
                break;
            case ':':
                current.ProcessLength();
                break;
            case '[':
                 current.ProcessComment();
                 break;
            default:
                // Anything unrecognized must be a tip name
                current.ProcessName();
                current.m_type = NodeType.TIP;
                current.m_corPt.y = (double)numTips + .5;
                treeTips.add(current);
                numTips++;
				}	
		}
		calcNodeCoords(root);
		//root.DumpTree(root);
		treeRead = true;
	}

	public void printTree() {
		PrinterJob pj = PrinterJob.getPrinterJob(); 
	    Book book = new Book();
	    PageFormat documentPageFormat = new PageFormat();
	    documentPageFormat.setOrientation(PageFormat.LANDSCAPE);
	    book.append(new Document(), documentPageFormat);
	    pj.setPageable(book);
	    if (pj.printDialog()) {
	        try {
	        	pj.print();
        	} catch (Exception PrintException) {
        		PrintException.printStackTrace();
	        }
	      }
	}
	
	private class Document implements Printable {
		public int print(Graphics g, PageFormat pageFormat, int page) {
			Graphics2D g2d = (Graphics2D) g;
			double scaleX = pageFormat.getImageableWidth()/frmRetreeControls.getWidth();
			double scaleY = pageFormat.getImageableHeight()/frmRetreeControls.getHeight();
			double plotScale = Math.min(scaleX, scaleY);
			g2d.scale(plotScale, plotScale);
			// ***JRMkludge***
			// There is apparently a margin around the page that cannot be drawn in. 
			// These offsets are to get out from under it. I suspect there is something
			// else I haven't allowed for in the graphics to page translation.
			// ***JRMkludge***
			int xoffset = 30;
			int yoffset = 50; 
			g2d.translate(pageFormat.getImageableX() + xoffset, pageFormat.getImageableY() + yoffset);
			Font titleFont = new Font("helvetica", Font.BOLD, 16);
			g2d.setFont(titleFont);
			g2d.drawString("Retree plot", xoffset, yoffset);
			plotRetree(g2d);
	      return (PAGE_EXISTS);
	    }
	}	
	
	public void saveTree() { 
	    JFileChooser fileChooser = new JFileChooser(filedir);
	    
	    // set up filter
	    FileNameExtensionFilter txtfilter = new FileNameExtensionFilter("Text file", "txt");
	    fileChooser.addChoosableFileFilter(txtfilter);
	    fileChooser.setSelectedFile(new File(outputTree));
	   
		int result = fileChooser.showSaveDialog(frmRetreeControls.getRootPane());
		if (result == JFileChooser.APPROVE_OPTION) {
			outputTree = fileChooser.getSelectedFile().getName();
			String curpath = fileChooser.getSelectedFile().getPath();
		    FileFilter ff = fileChooser.getFileFilter();
		    FileNameExtensionFilter extFilter = (FileNameExtensionFilter)ff;
		    String ext = extFilter.getExtensions()[0];
		    if (curpath.contains("." + ext))
		    {
			    outputTreePath = curpath;
		    }
		    else
		    {
			   outputTreePath = curpath + "." + ext;
		    }
			TestFileNames test = new TestFileNames();
			
			String opt = test.FileAlreadyExists(outputTreePath, "Tree file");
			if (opt != "q")
			{
				boolean append;
				if (opt == "w")
				{
					append = false;
				}
				else //if (opt == "a")
				{
					append = true;
				}

				String outputNewick = toNewick();
		        try {
		            BufferedWriter output = new BufferedWriter(new FileWriter(outputTreePath, append));
		            output.write(outputNewick);
		            output.close();
		        } catch ( IOException e ) {
		             e.printStackTrace();
		        }
			}
		}
	}
	
	public void showAbout() { 
		frmAbout = new JFrame();
		frmAbout.setVisible(true);
		frmAbout.setBounds(100, 50, 630, 300);
		frmAbout.setTitle("About Retree");
		frmAbout.setLayout(null);

		JLabel lblLine1 = new JLabel("Copyright 1993-2013. University of Washington and Joseph Felsenstein. All rights reserved.");
		lblLine1.setBounds(10, 10, 600, 14);
		frmAbout.add(lblLine1);
		
		JLabel lblLine2 = new JLabel("Permission is granted to reproduce, perform, and modify this program.");
		lblLine2.setBounds(10, 30, 600, 14);
		frmAbout.add(lblLine2);
   		
		JLabel lblLine3 = new JLabel("Permission is granted to distribute or provide access to this program provided that:");
		lblLine3.setBounds(10, 50, 600, 14);
		frmAbout.add(lblLine3);
   		
		JLabel lblLine4 = new JLabel("1) this copyright notice is not removed");
		lblLine4.setBounds(10,70, 600, 14);
		frmAbout.add(lblLine4);
   		
		JLabel lblLine5 = new JLabel("2) this program is not integrated with or called by any product or service that generates revenue");
		lblLine5.setBounds(10,90, 630, 14);
		frmAbout.add(lblLine5);
   		
		JLabel lblLine6 = new JLabel("3) your distribution of this program is free");
		lblLine6.setBounds(10,110, 600, 14);
		frmAbout.add(lblLine6);
   		
		JLabel lblLine7 = new JLabel("Any modified versions of this program that are distributed or accessible shall indicate");
		lblLine7.setBounds(10,150, 600, 14);
		frmAbout.add(lblLine7);
   		
		JLabel lblLine8 = new JLabel("that they are based on this program.  Educational institutions are granted permission");
		lblLine8.setBounds(10,170, 600, 14);
		frmAbout.add(lblLine8);
  		
		JLabel lblLine9 = new JLabel("to distribute this program to their students and staff for a fee to recover distribution costs.");
		lblLine9.setBounds(10,190, 600, 14);
		frmAbout.add(lblLine9);
  		
		JLabel lblLine10 = new JLabel("Permission requests for any other distribution of this program should be directed to:");
		lblLine10.setBounds(10,210, 600, 14);
		frmAbout.add(lblLine10);
  		
		JLabel lblLine11 = new JLabel("license (at) u.washington.edu.");
		lblLine11.setBounds(10,230, 600, 14);
		frmAbout.add(lblLine11);
		
		btnOK = new JButton("OK");
		btnOK.setBounds(500, 250, 84, 25);
		btnOK.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				frmAbout.dispose();
				frmRetreeControls.repaint();        					
			}
		});
		frmAbout.add(btnOK);
	}
	
	public void undoLast()
	{
		if (UndoTrees.size() > 1)
		{
			UndoTrees.remove(UndoTrees.size()-1);
			parseNewick(UndoTrees.get(UndoTrees.size()-1));
			calcNodeCoords(root);
			frmRetreeControls.repaint();
			if(UndoTrees.size() == 1)
			{
		        pntmUndo.setEnabled(false);						
		        mntmUndo.setEnabled(false);						
			}
		}		
	}
	
	public void saveInterumTree() { 
   		TestFileNames test = new TestFileNames();
		String opt = test.FileAlreadyExists(outputTreePath, "");
		boolean append = true;
		if (opt == "w")
		{
			append = false;
		}
		
		String outputNewick = toNewick();
        try {
            BufferedWriter output = new BufferedWriter(new FileWriter(outputTreePath,append));
            output.write(outputNewick);
            output.close();
        } catch ( IOException ioerr ) {
             ioerr.printStackTrace();
        }        
	}
	
	public String toNewick()
	{
		String retstr = "";
		retstr = addToNewick(root, retstr);
		
		return retstr;
	}
	
	public String addToNewick(NewickNode node, String retstr)
	{
		if(node.m_type == NodeType.TIP)
		{
			retstr += node.m_name;
			if (lengthFound)
			{
				retstr += ":"+node.m_length;
			}
		}
		else if (node.m_type == NodeType.ROOT)
		{
			Iterator<NewickNode> itr = node.m_children.iterator();
			retstr += "(";
			int numnode = 0;
			while(itr.hasNext())
			{
				if (numnode > 0)
				{
					retstr += ",";
				}
				retstr = addToNewick(itr.next(), retstr);
				numnode++;
			}
			retstr += ");\n";
		}
		else //node.m_type == NodeType.NODE
		{
			Iterator<NewickNode> itr = node.m_children.iterator();
			retstr += "(";
			int numnode = 0;
			while(itr.hasNext())
			{
				if (numnode > 0)
				{
					retstr += ",";
				}
				retstr = addToNewick(itr.next(), retstr);
				numnode++;
			}
			retstr += ")";
			if (lengthFound)
			{
				retstr += ":"+node.m_length;
			}
		}
		return retstr;
	}
	
	public double sumTreeLen (NewickNode node)
	{
		if (node.GetParent() == null)
		{
			return (0.0);
		}
		else 
		{
			return(sumTreeLen(node.GetParent()) + node.m_length);
		}
	}
	
	public void calcNodeCoords(NewickNode root)
	{
		// first fill in xcor
		walkDownFromRoot(root);
		
		// now do ycor
		for(int i=0; i<treeTips.size(); i++)
		{
			walkUpFromTip(treeTips.get(i));
		}
	}
	
	public void walkDownFromRoot(NewickNode node)
	{
		// calculates x coordinates
		node.m_corPt.x = 0.0;
		NewickNode parent = node.GetParent();
		if (parent != null)
		{
			node.m_corPt.x = parent.m_corPt.x + node.m_length;
		}

		if (node.m_type == NodeType.TIP)
		{
			return;
		}
		else
		{
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				walkDownFromRoot(itr.next());
			}
		}
	}
	
	public void walkUpFromTip(NewickNode tip)
	{
		// calculates y coordinates
		NewickNode parent = tip.GetParent();
		if (parent != null)
		{
			if (parent.m_corPt.y != 0.0)
			{
				walkUpFromTip(parent);
			}
			else
			{
				walkDownFromNode(parent);
			}
		}
	}
	
	public void walkDownFromNode(NewickNode node)
	{
		Iterator<NewickNode> itr = node.m_children.iterator();
		double ySum = 0.0;
		while(itr.hasNext())
		{
			double ynew = findLowestChildY(itr.next());
			ySum += ynew;
		}
		node.m_corPt.y = ySum/node.m_children.size();
	}
	
	public double findLowestChildY(NewickNode node)
	{
		if (node.m_type != NodeType.TIP)
		{
			if (node.m_corPt.y == 0.0)
			{
				walkDownFromNode(node);
			}
		}
		return node.m_corPt.y;
	}
	
	public void moveNodeBefore(NewickNode movingNode, NewickNode moveBeforeNode)
	{
		if ((movingNode.m_parent == moveBeforeNode.m_parent) && (movingNode.m_parent.m_children.size() == 2))
		{
			String msg = "Moving a node before a node that shares the moving node's parent is not meaningful.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);

			return;
		}
		if ((movingNode.m_parent == root) && (root.m_children.size() == 2))
		{	
			// root case with only 2 children is different because the root moves
			NewickNode newRoot = new NewickNode();
			Iterator<NewickNode> itr = root.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				if (child != movingNode)
				{
					movingNode.m_length += child.m_length;
					newRoot = child;
					
					// clear things out
					newRoot.m_length = 1.0;
					newRoot.m_type = NodeType.ROOT;
					newRoot.m_parent = null;
					newRoot.m_dispPt = new Point(0,0);
					newRoot.m_corPt = new Point2D.Double(0.0, 0.0);
				}
			}

			NewickNode newNode = moveBeforeNode.m_parent.AddChild();
			movingNode.m_parent = newNode;
			disconnectNode(moveBeforeNode);
			moveBeforeNode.m_parent = newNode;
			moveBeforeNode.m_length /= 2;
			newNode.m_length = moveBeforeNode.m_length;
			newNode.m_corPt.x = newNode.m_parent.m_corPt.x + newNode.m_length;
			newNode.m_children.add(moveBeforeNode);
			newNode.m_children.add(movingNode);
			root  = newRoot;
		}
		else
		{
			NewickNode newNode = moveBeforeNode.m_parent.AddChild();
			disconnectNode(movingNode);
			movingNode.m_parent = newNode;
			disconnectNode(moveBeforeNode);
			moveBeforeNode.m_parent = newNode;
			moveBeforeNode.m_length /= 2;
			newNode.m_length = moveBeforeNode.m_length;
			newNode.m_corPt.x = newNode.m_parent.m_corPt.x + newNode.m_length;
			newNode.m_children.add(moveBeforeNode);
			newNode.m_children.add(movingNode);
		}
		rebuildTipList();
		recalcYcors();
		calcNodeCoords(root);
	}
	
	public void moveNodeTo(NewickNode movingNode, NewickNode moveToNode)
	{
		if (movingNode.m_parent == moveToNode)
		{
			String msg = "Moving a node to it's own parent is not meaningful.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);

			return;
		}
		
		if (movingNode.m_parent == root) 
		{	
			// root case with only 2 children is different because the root moves
			NewickNode newRoot = new NewickNode();
			Iterator<NewickNode> itr = root.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				if (child != movingNode)
				{
					movingNode.m_length += child.m_length;
					newRoot = child;
					
					// clear things out
					newRoot.m_length = 1.0;
					newRoot.m_type = NodeType.ROOT;
					newRoot.m_parent = null;
					newRoot.m_dispPt = new Point(0,0);
					newRoot.m_corPt = new Point2D.Double(0.0, 0.0);
				}
			}
			disconnectNode(movingNode);
			movingNode.m_parent = moveToNode;
			moveToNode.m_children.add(movingNode);
			root  = newRoot;
		}
		else
		{
			disconnectNode(movingNode);
			movingNode.m_parent = moveToNode;
			moveToNode.m_children.add(movingNode);
		}
		rebuildTipList();
		recalcYcors();
		calcNodeCoords(root);
	}
	
	public void disconnectNode(NewickNode movingNode)
	{
		if (movingNode.m_parent.m_children.size() > 2)
		{
			// remove child
			Iterator<NewickNode> itr = movingNode.m_parent.m_children.iterator();
			while(itr.hasNext())
			{
				if (itr.next() == movingNode)
				{
					itr.remove();
				}
			}
		}
		else if (movingNode.m_parent.m_children.size() == 2)
		{
			if (movingNode.m_parent.m_type == NodeType.ROOT)
			{
				Iterator<NewickNode> itr = movingNode.m_parent.m_children.iterator();
				while(itr.hasNext())
				{
					NewickNode child = itr.next();
					if (child == movingNode)
					{
						itr.remove();
					}
				}		
			}
			else
			{
				// removing one child also eliminates the parent
				// transfer the length and the not moving child to the grandparent
				Iterator<NewickNode> itr = movingNode.m_parent.m_children.iterator();
				while(itr.hasNext())
				{
					NewickNode child = itr.next();
					if (child != movingNode)
					{
						// if .m_parent.m_parent == null, at the root and the rules change
						child.m_length += movingNode.m_parent.m_length;
						child.m_parent = movingNode.m_parent.m_parent;
						movingNode.m_parent.m_parent.m_children.add(child);
					}
				}	
				
				// eliminate the parent
				itr = movingNode.m_parent.m_parent.m_children.iterator();
				while(itr.hasNext())
				{
					NewickNode child = itr.next();
					if (child == movingNode.m_parent)
					{
						itr.remove();
						break;
					}
				}
			}
		}
		else
		{
			System.out.println("Somehow parent " + movingNode.m_parent.m_nodenum + " slipped through with a single child.");
		}		
	}

	public void transposeBranch(NewickNode node)
	{
		Collections.reverse(node.m_children);
		recalcYcors();
	}
	
	public void recalcYcors()
	{
		tipYcor = .5;
		resetYcors(root);
		for(int i=0; i<treeTips.size(); i++)
		{
			walkUpFromTip(treeTips.get(i));
		}		
	}
	
	public void resetYcors(NewickNode node)
	{
		// this prepares the tree so the internal nodes will be recalculated
		if (node.m_type == NodeType.TIP)
		{
			node.m_corPt.y = tipYcor;
			tipYcor += 1.0;
			return;
		}
		else
		{
			node.m_corPt.y = 0.0;
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				resetYcors(itr.next());
			}
		}		
	}
	
	public void rebuildTipList()
	{
		treeTips.clear();
		Iterator<NewickNode> itr = root.m_children.iterator();
		while(itr.hasNext())
		{
			NewickNode child = itr.next();
			findNextTip(child);
		}		
	}
	
	public void findNextTip(NewickNode node)
	{
		if (node.m_type == NodeType.TIP)
		{
			treeTips.add(node);
			return;
		}
		else
		{
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				findNextTip(itr.next());
			}		
		}
	}

	public void flipBranch(NewickNode node)
	{
		Collections.reverse(node.m_children);
		Iterator<NewickNode> itr = node.m_children.iterator();
		while(itr.hasNext())
		{
			NewickNode child = itr.next();
			flipBranch(child);
		}		
		recalcYcors();
	}
	
	public void addToUndo(String newTree)
	{
		UndoTrees.add(newTree);
		if (UndoTrees.size() > 1)
		{
			pntmUndo.setEnabled(true);
			mntmUndo.setEnabled(true);
		}
		else
		{
	        pntmUndo.setEnabled(false);
	        mntmUndo.setEnabled(false);
		}
	}
	
	public void outgroupRoot(NewickNode movingNode)
	{
		disconnectNode(movingNode);
		movingNode.m_length/=2;
		NewickNode newRoot = new NewickNode();
		//newRoot.NewRoot(newRoot);
		nodeNum++;
		newRoot.m_nodenum = nodeNum;
		newRoot.m_type = NodeType.ROOT;
		movingNode.m_parent = newRoot;
		
		NewickNode oldRoot = root;
		oldRoot.m_type = NodeType.NODE;
		oldRoot.m_parent = newRoot;
		oldRoot.m_length = movingNode.m_length;

		newRoot.m_children.add(oldRoot);
		newRoot.m_children.add(movingNode);
		root = newRoot;
		
		recalcYcors();
		calcNodeCoords(root);
	}
	
	public void midpointRoot()
	{
		// find the longest distance for each child of the root to its tip
		ArrayList<NewickNode> addedChildren= new ArrayList<NewickNode>();
		
		double delta = moveDist(root, addedChildren);
		
		// test for nodes inside delta
		Iterator<NewickNode> itr = root.m_children.iterator();
		newBase = root;
		boolean nodeInsideDelta = false;
		while(itr.hasNext())
		{
			NewickNode child = itr.next();
			if (child.m_length < delta)
			{
				nodeInsideDelta = true;
			}
		}

		if (nodeInsideDelta)
		{
			// decide which way to go
			itr = root.m_children.iterator();
			double maxdist = 0.0;
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				double dist = maxChildDist(child);
				if (dist > maxdist)
				{
					maxdist = dist;
					newBase = child;
				}
			}
		}

		if (newBase == root)
		{
			// simple case, adjust the lengths
			itr = root.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				if (child == longestBranch)
				{
					child.m_length -= delta;
				}
				else
				{
					child.m_length += delta;				
				}
			}
		}
		else
		{
			// unroot the tree, consolidating everything to newBase
			addedChildren.clear();
			itr = root.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				if (child != newBase)
				{
					child.m_length += newBase.m_length;
					child.m_parent = newBase;
					newBase.m_children.add(child);
					addedChildren.add(child);
				}
			}
			newBase.m_length = 0.0;
			newBase.m_parent = null;
			
			// recurse down the tree until no more nodes inside delta
			delta = moveDist(newBase, addedChildren);
			
			// now hook unrooted tree back to the root:			
			// remove longest branch from newBase
			itr = newBase.m_children.iterator();
			while(itr.hasNext())
			{
				if(itr.next() == longestBranch)
				{
					itr.remove();
				}
			}
		
			// adjust lengths 
			longestBranch.m_length -= delta;
			newBase.m_length = delta;
			
			// hook up children and root
			longestBranch.m_parent = root;
			newBase.m_parent = root;
			root.m_children.clear();
			root.m_children.add(longestBranch);
			root.m_children.add(newBase);
		}
		
		// reorder the branches
		if (root.m_children.size() > 2)
		{
			ArrayList<NewickNode> newChildList = new ArrayList<NewickNode>();
			newChildList.add(longestBranch);
			itr = root.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				if ((child != longestBranch) &&
					(child != nextLongBranch))
				{
					newChildList.add(child);
				}
			}		
			newChildList.add(nextLongBranch);
			root.m_children = newChildList;
		}
		
		recalcYcors();
		calcNodeCoords(root);
	}
	
	public double moveDist(NewickNode node, ArrayList<NewickNode> addedChildren)
	{
		double delta = calcMoveDist(node);
		
		if (node == root)
		{
			// special case for root node
			return delta;
		}
		else
		{
			// see if any nodes within delta of this node
			Iterator<NewickNode> itr = node.m_children.iterator();
			NewickNode closestNode = null;
			double closeLen = 1000000.0;
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				for(Iterator<NewickNode> i = addedChildren.iterator(); i.hasNext();)
				{
					if (child != i.next()) // keep things from oscillating
					{
						if (child.m_type == NodeType.NODE)
						{
							// base can only be a node
							if (child.m_length < delta)
							{
								if (child.m_length < closeLen)
								{
									closestNode = child;
									closeLen = child.m_length;
								}
							}
						}
					}
				}
			}
			
			if (closestNode == null)
			{
				// last recursion
				return calcMoveDist(newBase);
			}
			else
			{
				// shift base to closest node
				newBase = closestNode;
				addedChildren.clear();
				itr = node.m_children.iterator();
				while(itr.hasNext())
				{
					if (itr.next() == newBase)
					{
						itr.remove();
					}
				}
				node.m_parent = newBase;
				newBase.m_parent = null;
				node.m_length = newBase.m_length;
				newBase.m_length = 0.0;
				newBase.m_children.add(node);
				addedChildren.add(node);
				
				return moveDist(closestNode, addedChildren);
			}
		}
	}
	
	public Double calcMoveDist(NewickNode node)
	{
		Iterator<NewickNode> itr = node.m_children.iterator();
		double lenLB = 0.0;
		double lenNB = 0.0;
		longestBranch = null;
		nextLongBranch = null;
		while(itr.hasNext())
		{
			NewickNode child = itr.next();
			double brlen = maxChildDist(child);
			if (brlen > lenNB)
			{
				if(brlen > lenLB)
				{
					lenNB = lenLB;
					lenLB = brlen;
					nextLongBranch = longestBranch;
					longestBranch = child;
				}
				else
				{
					lenNB = brlen;					
					nextLongBranch = child;
				}
			}	
		}		
		return Math.abs((lenLB - lenNB)/2);	
	}
	
	public Double maxChildDist(NewickNode node)
	{
		if (node.m_type == NodeType.TIP){
			return node.m_length;
		}
		else
		{
			ArrayList<Double> tipDists = new ArrayList<Double>();	
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				tipDists.add(maxChildDist(child));
			}		
		
			double maxdist = 0.0;
			for(int i=0; i<tipDists.size(); i++)
			{
				if (tipDists.get(i) > maxdist)
				{
					maxdist = tipDists.get(i);
				}
			}
			return node.m_length + maxdist;	
		}
	}
	
		public double distToRoot(NewickNode node, double dist)
		{
		if (node.m_type == NodeType.ROOT)
		{
			return dist;
		}
		else
		{
			dist += node.m_length;
			return distToRoot(node.m_parent, dist);			
		}
	}
	
	public void registerShortcuts()
	{
		ActionMap amap = displayPanel.getActionMap();
		InputMap imap = displayPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		amap.put("ctlA", showAbout);
		imap.put(ctlA, "ctlA");
		amap.put("ctlP", printTree);
		imap.put(ctlP, "ctlP");
		amap.put("ctlQ", quitRun);
		imap.put(ctlQ, "ctlQ");
		amap.put("ctlR", readTree);
		imap.put(ctlR, "ctlR");
		amap.put("ctlS", saveTree);
		imap.put(ctlS, "ctlS");
		amap.put("ctlZ", doUndo);
		imap.put(ctlZ, "ctlZ");
	}
	
	public class PrintTree extends AbstractAction{
		public PrintTree(){
		}
		public void actionPerformed(ActionEvent e){
	        printTree();
			doMove = false;
		}
	};
	
	public class ReadTree extends AbstractAction{
		public ReadTree(){
		}
		public void actionPerformed(ActionEvent e){
    		UndoTrees.clear();
			readTree();
			frmRetreeControls.repaint();        					
	        mntmSave.setEnabled(true);
	        mntmPrint.setEnabled(true);
	        pntmName.setEnabled(true);
	        doMove = false;
		}
	};
	
	public class SaveTree extends AbstractAction{
		public SaveTree(){
		}
		public void actionPerformed(ActionEvent e){
			saveTree();
			doMove = false;
		}
	};
	
	public class QuitRun extends AbstractAction{
		public QuitRun(){
		}
		public void actionPerformed(ActionEvent e){
            System.exit(0);
		}
	};
	
	public class UndoLast extends AbstractAction{
		public UndoLast(){
		}
		public void actionPerformed(ActionEvent e){
			undoLast();
			doMove = false;
		}
	};
	
	public class ShowAbout extends AbstractAction{
		public ShowAbout(){
		}
		public void actionPerformed(ActionEvent e){
			showAbout();
		}
	};
	
	public void Redisplay(boolean doUndo)
	{
		if (doUndo)
		{
			addToUndo(toNewick());
		}
		root.ClearTree(root);
		frmRetreeControls.repaint();
	}

	public class NewickNode{
		Integer m_nodenum;
		NewickNode m_parent;
		ArrayList<NewickNode> m_children;
		String m_name;
		Double m_length;
		NodeType m_type;
		Point m_dispPt;
		Point2D.Double m_corPt;	
		Color m_color;
		
		public NewickNode()
		{
			m_nodenum = -1;
			m_parent = null;
		    m_children = new ArrayList<NewickNode>();			
			m_name = "";
			m_length = 1.0;
			m_type = NodeType.UNDEF;
			m_dispPt = new Point(0,0);
			m_corPt = new Point2D.Double(0.0, 0.0);
			m_color = Color.BLACK;
		}
		
		public NewickNode AddChild() {
		    newNode = new NewickNode();
		    m_children.add(newNode);
		    nodeNum++;
		    newNode.m_nodenum = nodeNum;
		    newNode.m_parent = this;
		    newNode.m_type = NodeType.NODE;
		    return newNode;
		}
		
		public NewickNode NewRoot(NewickNode node) {
			nodeNum++;
			node.m_nodenum = nodeNum;
			node.m_type = NodeType.ROOT;
		    return node;
		}

		public NewickNode GetParent(){
			return m_parent;
		} 
		
		public NewickNode Terminate()
		{
		    assert(m_children.size() > 1);
		    return this;
		} 

		public void ProcessLength()
		{
			++m_curr_char;  // skip the colon
			
			String lenval = "";
			int lenstr = 0;
            while ((m_newick[m_curr_char] != ',') &&
             	   (m_newick[m_curr_char] != ')'))
            {
            	lenval += m_newick[m_curr_char];
            	++m_curr_char;
            	lenstr++;
            }
            if (lenstr > 0)
            {
            	m_length = Double.valueOf(lenval);
            	lengthFound = true;
            }
		} 

		public void ProcessName()
		{
			int nchar = 0;
	        while ((m_newick[m_curr_char] != ':') &&
	        	   (m_newick[m_curr_char] != ',') &&
             	   (m_newick[m_curr_char] != ')'))
				{
	        		m_name += m_newick[m_curr_char];
					++m_curr_char;
					nchar++;
				}
	        if (nchar > maxNameLen)
	        {
	        	maxNameLen = nchar;
	        }
		} 	
	
		public void ProcessComment()
		// NB: We don't do anything with the contents of the comment yet
		{
	        while (m_newick[m_curr_char] != ']') 
			{
				++m_curr_char;
			}
		}
		
		public void DumpTree(NewickNode node)
		{
			// Debug dump of tree outward from a node to System.out
			System.out.println("Node number: " + node.m_nodenum);
			if (node.m_type == NodeType.ROOT)
			{
				System.out.println("Root node");
			}
			else if (node.m_type == NodeType.NODE)
			{
				System.out.println("Internal node");
			}
			else if (node.m_type == NodeType.TIP)
			{
				System.out.println("Tip node Name: " + node.m_name);
			}
			else
			{
				System.out.println("ERROR: Unknown node type");
			}
			System.out.println("Length: " + node.m_length);
			System.out.println("Xcor: " + node.m_corPt.x);
			System.out.println("Ycor: " + node.m_corPt.y);
			int n = 1;
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				System.out.println("     Child: "+ n );
				n++;
				DumpTree(itr.next());
			}
		}
		
		public void PlotTree(NewickNode node, Graphics g2d)
		{
			if(node.m_type == NodeType.TIP)
			{
				// print name
				node.m_dispPt.x = xOffset + rootLen + (int)(node.m_corPt.x * xScaleFactor);
				node.m_dispPt.y = (int)(node.m_corPt.y * yScaleFactor);
				int labelStartX = node.m_dispPt.x + nameOffsetX;
				int labelStartY = node.m_dispPt.y + nameTweakY;
				g2d.setColor(node.m_color);
				g2d.drawString(node.m_name,  labelStartX,  labelStartY);
				//System.out.println("Node Num: " + node.m_nodenum + " Type: " + node.m_type + " Xdisp: " + node.m_dispPt.x + " Ydisp: " + node.m_dispPt.y);
			
				return;
			}
			else if (node.m_type == NodeType.ROOT)
			{
				// draw root
				node.m_dispPt.x = xOffset + rootLen;
				node.m_dispPt.y = (int)(node.m_corPt.y * yScaleFactor);
				g2d.setColor(node.m_color);
				g2d.drawLine(xOffset, node.m_dispPt.y, node.m_dispPt.x,  node.m_dispPt.y); // root stub - arbitrary length
				
				// this is for the branch selector code
				newseg = new BranchSeg();
				newseg.m_startnode = node.m_parent;
				newseg.m_endnode = node;
				newseg.m_hbranch = new Line2D.Float(xOffset, node.m_dispPt.y, node.m_dispPt.x,  node.m_dispPt.y);
				branchSegs.add(newseg);				
				
			}
			
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				// draw branch
				NewickNode child = itr.next();
				child.m_dispPt.x = xOffset + rootLen + (int)(child.m_corPt.x * xScaleFactor);
				child.m_dispPt.y = (int)(child.m_corPt.y * yScaleFactor);
				g2d.setColor(child.m_color);
				g2d.drawLine(child.GetParent().m_dispPt.x, child.GetParent().m_dispPt.y, child.GetParent().m_dispPt.x,  child.m_dispPt.y); // y offset line
				g2d.drawLine(child.GetParent().m_dispPt.x, child.m_dispPt.y, child.m_dispPt.x,  child.m_dispPt.y); // x branch
				
				if (showLens)
				{
					Integer lenX = child.GetParent().m_dispPt.x + (child.m_dispPt.x - child.GetParent().m_dispPt.x)/2 - 10;
					g2d.drawString(Double.toString(child.m_length), lenX,  child.m_dispPt.y - 5);
				}

				g2d.setColor(Color.BLACK);
				//System.out.println("Node Num: " + node.m_nodenum + " Type: " + node.m_type + " Xdisp: " + node.m_dispPt.x + " Ydisp: " + node.m_dispPt.y);

				// this is for the branch selector code
				newseg = new BranchSeg();
				newseg.m_startnode = child.m_parent;
				newseg.m_endnode = child;
				newseg.m_vbranch = new Line2D.Float(child.GetParent().m_dispPt.x, child.GetParent().m_dispPt.y, child.GetParent().m_dispPt.x,  child.m_dispPt.y);
				newseg.m_hbranch = new Line2D.Float(child.GetParent().m_dispPt.x, child.m_dispPt.y, child.m_dispPt.x,  child.m_dispPt.y);
				branchSegs.add(newseg);	
				
				PlotTree(child, g2d);
			}
		}
		
		public void ClearTree(NewickNode node)
		{
			// reset color to black after edit
			if(node.m_type == NodeType.TIP)
			{
				node.m_color = Color.BLACK;
				return;
			}
			else if (node.m_type == NodeType.ROOT)
			{
				node.m_color = Color.BLACK;			
			}
			
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				node.m_color = Color.BLACK;
				ClearTree(child);
			}
		}
	}
	
	public class BranchSeg{
		Line2D.Float m_hbranch;
		Line2D.Float m_vbranch;
		NewickNode m_startnode;
		NewickNode m_endnode;
		public BranchSeg()
		{
			m_hbranch = null;
			m_vbranch = null;
			m_startnode = null;
			m_endnode = null;
		}
		
		public double distanceToBranch(Point2D.Float pick)
		{
			double a = pick.x - m_hbranch.x1;
			double b = pick.y - m_hbranch.y1;
			double c = m_hbranch.x2 - m_hbranch.x1;
			double d = m_hbranch.y2 - m_hbranch.y1;
			
			double dot_prod = a*c + b*d;
			double len_sq = c*c + d*d;
			double param = dot_prod / len_sq;
			
			double xx;
			double yy;
			
			if (param < 0 || (m_hbranch.x1 == m_hbranch.x2 && m_hbranch.y1 == m_hbranch.y2))
			{
				xx = m_hbranch.x1;
				yy = m_hbranch.y1;
			}
			else if (param > 1)
			{
				xx = m_hbranch.x2;
				yy = m_hbranch.y2;
			}
			else
			{
				xx = m_hbranch.x1 + param*c;
				yy = m_hbranch.y1 + param*d;
			}
			double dx = pick.x - xx;
			double dy = pick.y - yy;
			double dist = Math.sqrt(dx*dx + dy*dy);
			
			if (m_vbranch != null)
			{
				a = pick.x - m_vbranch.x1;
				b = pick.y - m_vbranch.y1;
				c = m_vbranch.x2 - m_vbranch.x1;
				d = m_vbranch.y2 - m_vbranch.y1;
				
				dot_prod = a*c + b*d;
				len_sq = c*c + d*d;
				param = dot_prod / len_sq;
				
				if (param < 0 || (m_vbranch.x1 == m_vbranch.x2 && m_vbranch.y1 == m_vbranch.y2))
				{
					xx = m_vbranch.x1;
					yy = m_vbranch.y1;
				}
				else if (param > 1)
				{
					xx = m_vbranch.x2;
					yy = m_vbranch.y2;
				}
				else
				{
					xx = m_vbranch.x1 + param*c;
					yy = m_vbranch.y1 + param*d;
				}
				dx = pick.x - xx;
				dy = pick.y - yy;
				if (Math.sqrt(dx*dx + dy*dy) < dist)
				{
					dist = Math.sqrt(dx*dx + dy*dy);
				}		
			}
			
			return dist;
		}
	}
}
