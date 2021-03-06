
// MFCApplication2Dlg.h : 头文件
//

#pragma once

#include <vector>

#define maxPicNum 5

using namespace std;


// CMFCApplication2Dlg 对话框
class CMFCApplication2Dlg : public CDialogEx
{
// 构造
public:
	CMFCApplication2Dlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_MFCAPPLICATION2_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()

	CPen pen;
	CClientDC *curDC;
	CPoint m_StartPoint;
	vector<vector<vector<CPoint>>> picLines;
	vector<vector<CPoint>> lines;
	vector<CPoint> currentLine;

	CStatic* pWnd;
	int curPicIndex;
	int picId[2][3][3]; // model type; distance: near mid far; method: ori, RS, our;
	int picIndices[5];

	void drawLines(vector<vector<CPoint>> &curLines);
	void sowPic(void);

	void readLines(void);
	void writeLines(void);

public:
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);



	afx_msg void OnBnClickedButtonUndo();
	afx_msg void OnBnClickedButtonNext();
	afx_msg void OnBnClickedButtonPre();
};
