
// MFCApplication2Dlg.cpp : 实现文件
//

#include "stdafx.h"
#include "MFCApplication2.h"
#include "MFCApplication2Dlg.h"
#include "afxdialogex.h"

#include <atlimage.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CMFCApplication2Dlg 对话框



CMFCApplication2Dlg::CMFCApplication2Dlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CMFCApplication2Dlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);

}

void CMFCApplication2Dlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CMFCApplication2Dlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON_NEXT, &CMFCApplication2Dlg::OnBnClickedButtonNext)
END_MESSAGE_MAP()


// CMFCApplication2Dlg 消息处理程序

BOOL CMFCApplication2Dlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO:  在此添加额外的初始化代码



	curPicIndex = -1;

	pWnd1 = (CStatic*)GetDlgItem(IDC_STATIC_PIC); // 得到 Picture Control 句柄


	pWnd2 = (CStatic*)GetDlgItem(IDC_STATIC_PIC2); // 得到 Picture Control 句柄

	GetDlgItem(IDC_STATIC_PIC2)->ModifyStyle(SS_TYPEMASK, SS_OWNERDRAW);
	
	CStatic *pDC = (CStatic*)GetDlgItem(IDC_STATIC);





	CRect rectCtrl;
	pWnd1->GetWindowRect(rectCtrl);
	ScreenToClient(rectCtrl);

	top = rectCtrl.top;
	w = rectCtrl.Width();
	h = rectCtrl.Height();

	// model type; distance: near mid far; method: ori, RS, our;

	picId[0][0][0] = IDB_PNG9;
	picId[0][0][1] = IDB_PNG11;
	picId[0][0][2] = IDB_PNG10;
	sizeI[0][0][0] = 572;
	sizeI[0][0][1] = 571;
	picId[0][1][0] = IDB_PNG6;
	picId[0][1][1] = IDB_PNG8;
	picId[0][1][2] = IDB_PNG7;
	sizeI[0][1][0] = 286;
	sizeI[0][1][1] = 286;
	picId[0][2][0] = IDB_PNG3;
	picId[0][2][1] = IDB_PNG5;
	picId[0][2][2] = IDB_PNG4;
	sizeI[0][2][0] = 142;
	sizeI[0][2][1] = 142;


	picId[1][0][0] = IDB_PNG22;
	picId[1][0][1] = IDB_PNG24;
	picId[1][0][2] = IDB_PNG23;
	sizeI[1][0][0] = 349;
	sizeI[1][0][1] = 427;
	picId[1][1][0] = IDB_PNG19;
	picId[1][1][1] = IDB_PNG21;
	picId[1][1][2] = IDB_PNG20;
	sizeI[1][1][0] = 171;
	sizeI[1][1][1] = 207;
	picId[1][2][0] = IDB_PNG16;
	picId[1][2][1] = IDB_PNG18;
	picId[1][2][2] = IDB_PNG17;
	sizeI[1][2][0] = 85;
	sizeI[1][2][1] = 103;

	picId[2][0][0] = IDB_PNG1;
	picId[2][0][1] = IDB_PNG12;
	picId[2][0][2] = IDB_PNG2;
	sizeI[2][0][0] = 174;
	sizeI[2][0][1] = 290;
	picId[2][1][0] = IDB_PNG13;
	picId[2][1][1] = IDB_PNG15;
	picId[2][1][2] = IDB_PNG14;
	sizeI[2][1][0] = 86;
	sizeI[2][1][1] = 144;

	picId[3][0][0] = IDB_PNG25;
	picId[3][0][1] = IDB_PNG27;
	picId[3][0][2] = IDB_PNG26;
	sizeI[3][0][0] = 140;
	sizeI[3][0][1] = 143;
	picId[3][1][0] = IDB_PNG28;
	picId[3][1][1] = IDB_PNG30;
	picId[3][1][2] = IDB_PNG29;
	sizeI[3][1][0] = 279;
	sizeI[3][1][1] = 286;

	bool hasSelect[10] = { false };

	srand((unsigned int)(time(NULL)));

	for (int i = 0; i < 10; i++)
	{

		int picIndex = (int)((float)(rand() - 1) / RAND_MAX * 10);
		while (hasSelect[picIndex])
			picIndex = (int)((float)(rand() - 1) / RAND_MAX * 10);
		hasSelect[picIndex] = true;

		int mIndex;
		int dIndex;


		if (picIndex < 6)
		{
			mIndex = picIndex / 3;
			dIndex = picIndex % 3;
		}
		else
		{
			mIndex = picIndex / 2-1;
			dIndex = picIndex % 2;

		}


		bool meHas[3] = { false };
		for (int j = 0; j < 3; j++)
		{
			int msIndex = (int)((float)(rand() - 1) / RAND_MAX * 3);
			while (meHas[msIndex])
				msIndex = (int)((float)(rand() - 1) / RAND_MAX * 3);
			meHas[msIndex] = true;

			pos[i + j * 10][2] = sizeI[mIndex][dIndex][0];
			pos[i + j * 10][3] = sizeI[mIndex][dIndex][1];

			if ((int)((float)(rand() - 1) / RAND_MAX * 2) == 0)
			{
				picIndices[i + j * 10][0] = picId[mIndex][dIndex][msIndex];
				picIndices[i + j * 10][1] = picId[mIndex][dIndex][(msIndex + 1) % 3];
			}
			else
			{
				picIndices[i + j * 10][1] = picId[mIndex][dIndex][msIndex];
				picIndices[i + j * 10][0] = picId[mIndex][dIndex][(msIndex + 1) % 3];
			}
		}
	}




		  

	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void CMFCApplication2Dlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CMFCApplication2Dlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CMFCApplication2Dlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}




void CMFCApplication2Dlg::OnBnClickedButtonNext()
{
	// TODO:  在此添加控件通知处理程序代码

	int i = ((CButton *)GetDlgItem(IDC_LEFT))->GetCheck();
	int j = ((CButton *)GetDlgItem(IDC_RIGHT))->GetCheck();


	int ii = ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck();
	int jj = ((CButton *)GetDlgItem(IDC_RADIO2))->GetCheck();


	if (((i == 0 && j == 0) || (ii == 0 && jj == 0)) && curPicIndex >= 0)
		return;

	curPicIndex++;

	start = finish;
	finish = clock();

	if (curPicIndex - 1 >= 0)
	{

		timeLast[curPicIndex-1] = (double)(finish - start) / CLOCKS_PER_SEC;



		if (i == 1)
			result[curPicIndex - 1] = 0;
		else
			result[curPicIndex - 1] = 1;


		if (ii == 1)
			result2[curPicIndex - 1] = 0;
		else
			result2[curPicIndex - 1] = 1;

	}


	if (curPicIndex >= maxPicNum)
	{
		FILE *fout;
		fopen_s(&fout, "result.dat", "wt");
		while (!fout)
		{
			fopen_s(&fout, "result.dat", "wt");
		}
		fwrite(&picIndices[0], sizeof(int), 60, fout);

		fwrite(&result[0], sizeof(int), 30, fout);
		fwrite(&result2[0], sizeof(int), 30, fout);
		fwrite(&timeLast[0], sizeof(double), 30, fout);

		fclose(fout);
		exit(0);
	}
	sowPic();





	((CButton *)GetDlgItem(IDC_LEFT))->SetCheck(FALSE);//不选上
	((CButton *)GetDlgItem(IDC_RIGHT))->SetCheck(FALSE);//不选上

	((CButton *)GetDlgItem(IDC_RADIO1))->SetCheck(FALSE);//不选上
	((CButton *)GetDlgItem(IDC_RADIO2))->SetCheck(FALSE);//不选上
	
}

void CMFCApplication2Dlg::sowPic(void)
{


	//设置静态控件窗口风格为位图居中显示  
	pWnd1->ModifyStyle(0xf, SS_BITMAP|SS_CENTERIMAGE);
	pWnd1->SetWindowPos(NULL, 603 - pos[curPicIndex][2], 300 - pos[curPicIndex][3]/2, pos[curPicIndex][2], pos[curPicIndex][3], SWP_NOZORDER);

	CPngImage image1;
	image1.Load(picIndices[curPicIndex][0]);
	//显示图片
	pWnd1->SetBitmap((HBITMAP)image1);


	//设置静态控件窗口风格为位图居中显示  
	pWnd2->ModifyStyle(0xf, SS_BITMAP | SS_CENTERIMAGE);
	pWnd2->SetWindowPos(NULL, 620, 300 - pos[curPicIndex][3] / 2, pos[curPicIndex][2], pos[curPicIndex][3], SWP_NOZORDER);
	CPngImage image2;
	image2.Load(picIndices[curPicIndex][1]);
	//显示图片
	pWnd2->SetBitmap((HBITMAP)image2);
	


}


void CMFCApplication2Dlg::writeResult(void)
{

}
