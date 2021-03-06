
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
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_BN_CLICKED(IDC_BUTTON_UNDO, &CMFCApplication2Dlg::OnBnClickedButtonUndo)
	ON_BN_CLICKED(IDC_BUTTON_NEXT, &CMFCApplication2Dlg::OnBnClickedButtonNext)
	ON_BN_CLICKED(IDC_BUTTON_PRE, &CMFCApplication2Dlg::OnBnClickedButtonPre)
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



	COLORREF m_Color(RGB(255, 0, 0));
	pen.CreatePen(PS_SOLID, 5, m_Color);
	curDC = new CClientDC(this);
	curDC->SelectObject(&pen);
	curPicIndex = -1;
	picLines.resize(maxPicNum);



	pWnd = (CStatic*)GetDlgItem(IDC_STATIC_PIC); // 得到 Picture Control 句柄

	//设置静态控件窗口风格为位图居中显示  
	pWnd->ModifyStyle(0xf, SS_BITMAP | SS_CENTERIMAGE);

	// model type; distance: near mid far; method: ori, RS, our;

	picId[0][0][0] = IDB_PNG1;
	picId[0][0][1] = IDB_PNG12;
	picId[0][0][2] = IDB_PNG2;

	picId[0][1][0] = IDB_PNG13;
	picId[0][1][1] = IDB_PNG15;
	picId[0][1][2] = IDB_PNG14;


	picId[1][0][0] = IDB_PNG9;
	picId[1][0][1] = IDB_PNG11;
	picId[1][0][2] = IDB_PNG10;

	picId[1][1][0] = IDB_PNG6;
	picId[1][1][1] = IDB_PNG8;
	picId[1][1][2] = IDB_PNG7;

	picId[1][2][0] = IDB_PNG3;
	picId[1][2][1] = IDB_PNG5;
	picId[1][2][2] = IDB_PNG4;



		  

	srand((unsigned int)(time(NULL)));
	bool hasSelect[3] = { false };

	
	for (int index = 0; index < 2; index++)
	{
		int methIndex = (int)((float)(rand() - 1) / RAND_MAX * 3);

		while (hasSelect[methIndex])
		{
			methIndex = (int)((float)(rand() - 1) / RAND_MAX * 3);
		}
		hasSelect[methIndex] = true;

		picIndices[index] = picId[0][index][methIndex];
	}

	hasSelect[0] = false;
	hasSelect[1] = false;
	hasSelect[2] = false;

	for (int index = 0; index < 3; index++)
	{
		int methIndex = (int)((float)(rand() - 1) / RAND_MAX * 3);
		
		while (hasSelect[methIndex])
		{
			methIndex = (int)((float)(rand() - 1) / RAND_MAX * 3);
		}
		hasSelect[methIndex] = true;

		picIndices[index+2] = picId[1][index][methIndex];
	}
	

	readLines();


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


void CMFCApplication2Dlg::OnMouseMove(UINT nFlags, CPoint point)
{
	// TODO:  在此添加消息处理程序代码和/或调用默认值

	if (curPicIndex < 0)
		return;

	if (nFlags == MK_LBUTTON) //判断鼠标左键是否按下，如果按下，则移动时画线
	{
		curDC->MoveTo(m_StartPoint);
		curDC->LineTo(point);
		m_StartPoint = point; //将画线的起点移动到鼠标移动后的点

		currentLine.push_back(m_StartPoint);
	}


	CDialogEx::OnMouseMove(nFlags, point);
}


void CMFCApplication2Dlg::OnLButtonDown(UINT nFlags, CPoint point)
{
	// TODO:  在此添加消息处理程序代码和/或调用默认值

	if (curPicIndex < 0)
		return;
	m_StartPoint = point;
	currentLine.push_back(m_StartPoint);
	CDialogEx::OnLButtonDown(nFlags, point);
}


void CMFCApplication2Dlg::OnLButtonUp(UINT nFlags, CPoint point)
{
	// TODO:  在此添加消息处理程序代码和/或调用默认值

	if (curPicIndex < 0)
		return;
	curDC->MoveTo(m_StartPoint);
	curDC->LineTo(point);
	currentLine.push_back(point);
	picLines[curPicIndex].push_back(currentLine);
	currentLine.clear();

	CDialogEx::OnLButtonUp(nFlags, point);
}


void CMFCApplication2Dlg::OnBnClickedButtonUndo()
{
	// TODO:  在此添加控件通知处理程序代码

	RedrawWindow();
	if (picLines[curPicIndex].empty())
		return;
	picLines[curPicIndex].pop_back();
	drawLines(picLines[curPicIndex]);
}


void CMFCApplication2Dlg::OnBnClickedButtonNext()
{
	// TODO:  在此添加控件通知处理程序代码
	curPicIndex++;
	if (curPicIndex >= maxPicNum)
		curPicIndex--;
	writeLines();
	sowPic();
}

void CMFCApplication2Dlg::OnBnClickedButtonPre()
{
	// TODO:  在此添加控件通知处理程序代码

	curPicIndex--;
	if (curPicIndex < 0)
		curPicIndex++;
		
	writeLines();
	sowPic();


}


void CMFCApplication2Dlg::sowPic(void)
{


	CPngImage image;
	image.Load(picIndices[curPicIndex]);
	//显示图片
	pWnd->SetBitmap((HBITMAP)image);

	drawLines(picLines[curPicIndex]);

}

void CMFCApplication2Dlg::drawLines(vector<vector<CPoint>> &curLines)
{
	for (unsigned int lIndex = 0; lIndex < curLines.size(); lIndex++)
	{
		vector<CPoint> &curLin = curLines[lIndex];
		curDC->MoveTo(curLin[0]);
		for (unsigned int pIndex = 1; pIndex < curLin.size(); pIndex++)
		{
			curDC->LineTo(curLin[pIndex]);
		}
	}
}


void CMFCApplication2Dlg::readLines(void)
{
	FILE *fin;
	fopen_s(&fin, "result.dat", "rb");

	if (!fin)
		return;
	fread(&picIndices[0], sizeof(int), 3, fin);
	int pNum;
	fread(&pNum, sizeof(int), 1, fin);
	picLines.resize(pNum);
	for (int pIndex = 0; pIndex < pNum; pIndex++)
	{
		int lNum;
		fread(&lNum, sizeof(int), 1, fin);
		picLines[pIndex].resize(lNum);
		for (int lIndex = 0; lIndex < lNum; lIndex++)
		{
			int cNum;
			fread(&cNum, sizeof(int), 1, fin);
			picLines[pIndex][lIndex].resize(cNum);
			fread(&picLines[pIndex][lIndex][0], sizeof(int), cNum*2, fin);
		}
	}

	fclose(fin);
}
void CMFCApplication2Dlg::writeLines(void)
{

	FILE *fout;
	fopen_s(&fout, "result.dat", "wb");
	while (!fout)
	{
		fopen_s(&fout, "result.dat", "wb");
	}
	fwrite(&picIndices[0], sizeof(int), 3, fout);
	int pNum = picLines.size();
	fwrite(&pNum, sizeof(int), 1, fout);
	for (int pIndex = 0; pIndex < pNum; pIndex++)
	{
		int lNum = picLines[pIndex].size();
		fwrite(&lNum, sizeof(int), 1, fout);
		for (int lIndex = 0; lIndex < lNum; lIndex++)
		{
			int cNum = picLines[pIndex][lIndex].size();
			fwrite(&cNum, sizeof(int), 1, fout);
			fwrite(&picLines[pIndex][lIndex][0], sizeof(int), cNum*2, fout);
		}
	}
	fclose(fout);
}
