
// MFCApplication2Dlg.h : 头文件
//

#pragma once

#include <vector>

#define maxPicNum 15

#include "time.h"

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
	

	void writeResult(void);

	CStatic* pWnd1;
	CStatic* pWnd2;
	int curPicIndex;
	int picId[4][3][3]; // model type; distance: near mid far; method: ori, RS, our;
	int picIndices[30][2];
	int pos[30][4];
	int sizeI[4][3][2];
	int w, h, top;

	int result[31];
	int result2[31];
	double timeLast[31];
	void sowPic(void);

	clock_t start, finish;

public:


	afx_msg void OnBnClickedButtonNext();
};
