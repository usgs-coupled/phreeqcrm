// SelectedOutput.cpp: implementation of the CSelectedOutput class.
//
//////////////////////////////////////////////////////////////////////
#if defined(WIN32)
#include <windows.h>            // reqd to avoid namespace problems
#endif

#if defined(_DEBUG)
#include <sstream>              // std::ostringstream
#endif

#include <stdarg.h>
#include <stdio.h>

#if defined(PHREEQC_CLASS)
#include "phrqtype.h"
#include "p2c.h"
#include "global_structures.h"
#include "basic.h"
#include "Phreeqc.h"

// COMMENT: {2/24/2010 6:01:56 PM}extern int user_punch_count_headings;
// COMMENT: {2/24/2010 6:01:56 PM}extern char **user_punch_headings;
#endif

#include "SelectedOutput.hxx"
#include "phreeqcns.hxx"

const size_t RESERVE_ROWS = 80;
const size_t RESERVE_COLS = 80;

// COMMENT: {3/3/2010 5:31:34 PM}int EndRow(void);
void AddSelectedOutput(const char* name, const char* format, va_list argptr);
int warning_msg (const char *err_str);

// COMMENT: {3/3/2010 8:55:29 PM}int Phreeqc::EndRow(void)
// COMMENT: {3/3/2010 8:55:29 PM}{
// COMMENT: {3/3/2010 8:55:29 PM}// COMMENT: {3/3/2010 7:29:42 PM}	if (CSelectedOutput::Instance()->GetRowCount() <= 1)
// COMMENT: {3/3/2010 8:55:29 PM}	if (this->SelectedOutput.GetRowCount() <= 1)
// COMMENT: {3/3/2010 8:55:29 PM}	{
// COMMENT: {3/3/2010 8:55:29 PM}		// ensure all user_punch headings are included
// COMMENT: {3/3/2010 8:55:29 PM}		for (int i = n_user_punch_index; i < user_punch_count_headings; ++i)
// COMMENT: {3/3/2010 8:55:29 PM}		{
// COMMENT: {3/3/2010 8:55:29 PM}			CSelectedOutput::Instance()->PushBackEmpty(user_punch_headings[i]);
// COMMENT: {3/3/2010 8:55:29 PM}		}
// COMMENT: {3/3/2010 8:55:29 PM}	}
// COMMENT: {3/3/2010 8:55:29 PM}	return CSelectedOutput::Instance()->EndRow();
// COMMENT: {3/3/2010 8:55:29 PM}}

// COMMENT: {3/3/2010 8:55:40 PM}int PushBackDouble(const char* key, double dVal)
// COMMENT: {3/3/2010 8:55:40 PM}{
// COMMENT: {3/3/2010 8:55:40 PM}	return CSelectedOutput::Instance()->PushBackDouble(key, dVal);
// COMMENT: {3/3/2010 8:55:40 PM}}
// COMMENT: {3/3/2010 8:55:40 PM}
// COMMENT: {3/3/2010 8:55:40 PM}int PushBackLong(const char* key, long lVal)
// COMMENT: {3/3/2010 8:55:40 PM}{
// COMMENT: {3/3/2010 8:55:40 PM}	return CSelectedOutput::Instance()->PushBackLong(key, lVal);
// COMMENT: {3/3/2010 8:55:40 PM}}
// COMMENT: {3/3/2010 8:55:40 PM}
// COMMENT: {3/3/2010 8:55:40 PM}int PushBackString(const char* key, const char* sVal)
// COMMENT: {3/3/2010 8:55:40 PM}{
// COMMENT: {3/3/2010 8:55:40 PM}	return CSelectedOutput::Instance()->PushBackString(key, sVal);
// COMMENT: {3/3/2010 8:55:40 PM}}
// COMMENT: {3/3/2010 8:55:40 PM}
// COMMENT: {3/3/2010 8:55:40 PM}int PushBackEmpty(const char* key)
// COMMENT: {3/3/2010 8:55:40 PM}{
// COMMENT: {3/3/2010 8:55:40 PM}	return CSelectedOutput::Instance()->PushBackEmpty(key);
// COMMENT: {3/3/2010 8:55:40 PM}}


// COMMENT: {3/3/2010 8:58:25 PM}void AddSelectedOutput(const char* name, const char* format, va_list argptr)
// COMMENT: {3/3/2010 8:58:25 PM}{
// COMMENT: {3/3/2010 8:58:25 PM}	int bInt;
// COMMENT: {3/3/2010 8:58:25 PM}	int bDouble;
// COMMENT: {3/3/2010 8:58:25 PM}	int bString;
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}	int state;
// COMMENT: {3/3/2010 8:58:25 PM}	int bLongDouble;
// COMMENT: {3/3/2010 8:58:25 PM}	char ch;
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}	/* state values
// COMMENT: {3/3/2010 8:58:25 PM}	0 Haven't found start(%)
// COMMENT: {3/3/2010 8:58:25 PM}	1 Just read start(%)
// COMMENT: {3/3/2010 8:58:25 PM}	2 Just read Flags(-0+ #) (zero or more)
// COMMENT: {3/3/2010 8:58:25 PM}	3 Just read Width
// COMMENT: {3/3/2010 8:58:25 PM}	4 Just read Precision start (.)
// COMMENT: {3/3/2010 8:58:25 PM}	5 Just read Size modifier
// COMMENT: {3/3/2010 8:58:25 PM}	6 Just read Type
// COMMENT: {3/3/2010 8:58:25 PM}	*/
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}	if (name == NULL) {
// COMMENT: {3/3/2010 8:58:25 PM}		return;
// COMMENT: {3/3/2010 8:58:25 PM}	}
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}	bDouble = 0;
// COMMENT: {3/3/2010 8:58:25 PM}	bInt = 0;
// COMMENT: {3/3/2010 8:58:25 PM}	bString = 0;
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}	bLongDouble = 0;
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}	state = 0;
// COMMENT: {3/3/2010 8:58:25 PM}	ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}	while (ch != '\0') {
// COMMENT: {3/3/2010 8:58:25 PM}		switch (state) {
// COMMENT: {3/3/2010 8:58:25 PM}	case 0: /* looking for Start specification (%) */
// COMMENT: {3/3/2010 8:58:25 PM}		switch (ch) {
// COMMENT: {3/3/2010 8:58:25 PM}	case '%':
// COMMENT: {3/3/2010 8:58:25 PM}		state = 1;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	default:
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}		ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 1: /* reading Flags (zero or more(-,+,0,# or space)) */
// COMMENT: {3/3/2010 8:58:25 PM}		switch (ch) {
// COMMENT: {3/3/2010 8:58:25 PM}	case '-': case '0': case '+': case ' ': case '#':
// COMMENT: {3/3/2010 8:58:25 PM}		ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	default:
// COMMENT: {3/3/2010 8:58:25 PM}		state = 2;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 2: /* reading Minimum field width (decimal integer constant) */
// COMMENT: {3/3/2010 8:58:25 PM}		switch (ch) {
// COMMENT: {3/3/2010 8:58:25 PM}	case '.':
// COMMENT: {3/3/2010 8:58:25 PM}		state = 3;
// COMMENT: {3/3/2010 8:58:25 PM}		ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
// COMMENT: {3/3/2010 8:58:25 PM}		ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	default:
// COMMENT: {3/3/2010 8:58:25 PM}		state = 4;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 3: /* reading Precision specification (period already read) */
// COMMENT: {3/3/2010 8:58:25 PM}		switch (ch)
// COMMENT: {3/3/2010 8:58:25 PM}		{
// COMMENT: {3/3/2010 8:58:25 PM}		case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
// COMMENT: {3/3/2010 8:58:25 PM}			ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}			break;
// COMMENT: {3/3/2010 8:58:25 PM}		default:
// COMMENT: {3/3/2010 8:58:25 PM}			state = 4;
// COMMENT: {3/3/2010 8:58:25 PM}			break;
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 4: /* reading Size modifier */
// COMMENT: {3/3/2010 8:58:25 PM}		switch (ch)
// COMMENT: {3/3/2010 8:58:25 PM}		{
// COMMENT: {3/3/2010 8:58:25 PM}		case 'l':
// COMMENT: {3/3/2010 8:58:25 PM}			ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}			break;
// COMMENT: {3/3/2010 8:58:25 PM}		case 'L':
// COMMENT: {3/3/2010 8:58:25 PM}			bLongDouble = 1;
// COMMENT: {3/3/2010 8:58:25 PM}			ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}			break;
// COMMENT: {3/3/2010 8:58:25 PM}		case 'h':
// COMMENT: {3/3/2010 8:58:25 PM}			ch = *format++;
// COMMENT: {3/3/2010 8:58:25 PM}			break;
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}		state = 5;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 5: /* reading Conversion letter */
// COMMENT: {3/3/2010 8:58:25 PM}		switch (ch) {
// COMMENT: {3/3/2010 8:58:25 PM}	case 'c':
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 'd':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'i':
// COMMENT: {3/3/2010 8:58:25 PM}		bInt = 1;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 'n':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'o':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'p':
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 's':
// COMMENT: {3/3/2010 8:58:25 PM}		bString = 1;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 'u':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'x':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'X':
// COMMENT: {3/3/2010 8:58:25 PM}	case '%':
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	case 'f':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'e':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'E':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'g':
// COMMENT: {3/3/2010 8:58:25 PM}	case 'G':
// COMMENT: {3/3/2010 8:58:25 PM}		bDouble = 1;
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}	default:
// COMMENT: {3/3/2010 8:58:25 PM}		ASSERT(false);
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}		ch = '\0';  /* done */
// COMMENT: {3/3/2010 8:58:25 PM}		break;
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}	}
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}	if (bDouble) {
// COMMENT: {3/3/2010 8:58:25 PM}		double valDouble;
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}		if (bLongDouble) {
// COMMENT: {3/3/2010 8:58:25 PM}			valDouble = (double)va_arg(argptr, long double);
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}		else {
// COMMENT: {3/3/2010 8:58:25 PM}			valDouble = va_arg(argptr, double);
// COMMENT: {3/3/2010 8:58:25 PM}		}
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}		CSelectedOutput::Instance()->PushBackDouble(name, valDouble);
// COMMENT: {3/3/2010 8:58:25 PM}	}
// COMMENT: {3/3/2010 8:58:25 PM}	else if (bInt) {
// COMMENT: {3/3/2010 8:58:25 PM}		int valInt;
// COMMENT: {3/3/2010 8:58:25 PM}		valInt = va_arg(argptr, int);
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}		CSelectedOutput::Instance()->PushBackLong(name, (long)valInt);
// COMMENT: {3/3/2010 8:58:25 PM}	}
// COMMENT: {3/3/2010 8:58:25 PM}	else if (bString) {
// COMMENT: {3/3/2010 8:58:25 PM}		char* valString;
// COMMENT: {3/3/2010 8:58:25 PM}		valString = (char *)va_arg(argptr, char *);
// COMMENT: {3/3/2010 8:58:25 PM}
// COMMENT: {3/3/2010 8:58:25 PM}		CSelectedOutput::Instance()->PushBackString(name, valString);
// COMMENT: {3/3/2010 8:58:25 PM}	}
// COMMENT: {3/3/2010 8:58:25 PM}	else {
// COMMENT: {3/3/2010 8:58:25 PM}		ASSERT(false);
// COMMENT: {3/3/2010 8:58:25 PM}		CSelectedOutput::Instance()->PushBackEmpty(name);
// COMMENT: {3/3/2010 8:58:25 PM}	}
// COMMENT: {3/3/2010 8:58:25 PM}}

// COMMENT: {3/3/2010 8:56:03 PM}// COMMENT: {11/16/2004 10:18:22 PM}CSelectedOutput CSelectedOutput::singleton;
// COMMENT: {3/3/2010 8:56:03 PM}CSelectedOutput* CSelectedOutput::s_instance = 0;
// COMMENT: {3/3/2010 8:56:03 PM}
// COMMENT: {3/3/2010 8:56:03 PM}CSelectedOutput* CSelectedOutput::Instance()
// COMMENT: {3/3/2010 8:56:03 PM}{
// COMMENT: {3/3/2010 8:56:03 PM}	if (s_instance == 0)
// COMMENT: {3/3/2010 8:56:03 PM}	{
// COMMENT: {3/3/2010 8:56:03 PM}		s_instance = new CSelectedOutput;
// COMMENT: {3/3/2010 8:56:03 PM}	}
// COMMENT: {3/3/2010 8:56:03 PM}	return s_instance;
// COMMENT: {3/3/2010 8:56:03 PM}}
// COMMENT: {3/3/2010 8:56:03 PM}
// COMMENT: {3/3/2010 8:56:03 PM}void CSelectedOutput::Release()
// COMMENT: {3/3/2010 8:56:03 PM}{
// COMMENT: {3/3/2010 8:56:03 PM}	if (s_instance)
// COMMENT: {3/3/2010 8:56:03 PM}	{
// COMMENT: {3/3/2010 8:56:03 PM}		delete s_instance;
// COMMENT: {3/3/2010 8:56:03 PM}		s_instance = 0;
// COMMENT: {3/3/2010 8:56:03 PM}	}
// COMMENT: {3/3/2010 8:56:03 PM}}

CSelectedOutput::CSelectedOutput()
: m_nRowCount(0)
{
	this->m_arrayVar.reserve(RESERVE_COLS);
}

CSelectedOutput::~CSelectedOutput()
{
}

void CSelectedOutput::Clear(void)
{
	this->m_nRowCount = 0;
	this->m_vecVarHeadings.clear();
	this->m_arrayVar.clear();
	this->m_mapHeadingToCol.clear();
}

size_t CSelectedOutput::GetRowCount(void)const
{
	if (this->GetColCount())
	{
		return this->m_nRowCount + 1; // rows + heading
	}
	else
	{
		return 0;
	}
}

size_t CSelectedOutput::GetColCount(void)const
{
	return this->m_vecVarHeadings.size();
}

CVar CSelectedOutput::Get(int nRow, int nCol)const
{
	CVar v;
	this->Get(nRow, nCol, &v);
	return v;
}

VRESULT CSelectedOutput::Get(int nRow, int nCol, VAR* pVAR)const
{
	if ((size_t)nRow >= this->GetRowCount() || (size_t)nRow < 0) {
		pVAR->type = TT_ERROR;
		pVAR->vresult = VR_INVALIDROW;
		return pVAR->vresult;
	}
	if ((size_t)nCol >= this->GetColCount() || (size_t)nCol < 0) {
		pVAR->type = TT_ERROR;
		pVAR->vresult = VR_INVALIDCOL;
		return pVAR->vresult;
	}
	if (nRow)
	{
		ASSERT((size_t)nRow <= this->m_arrayVar[nCol].size());
		return ::VarCopy(pVAR, &(this->m_arrayVar[nCol])[nRow - 1]);
	}
	else
	{
		return ::VarCopy(pVAR, &(this->m_vecVarHeadings[nCol]));
	}
}


int CSelectedOutput::EndRow(void)
{
	if (size_t ncols = this->GetColCount())
	{
		++this->m_nRowCount;

		// make sure array is full
		for (size_t col = 0; col < ncols; ++col)
		{
			for (size_t row = 0; row < this->m_nRowCount; ++row)
			{
				size_t nrows = this->m_arrayVar[col].size();
				if (nrows < this->m_nRowCount)
				{
					// fill w/ empty
					this->m_arrayVar[col].resize(this->m_nRowCount);
				}
#if defined(_DEBUG)
				else if (nrows > this->m_nRowCount)
				{
					ASSERT(false);
				}
#endif
			}
		}
	}
// COMMENT: {11/27/2006 7:22:37 PM}#if defined(_DEBUG)
// COMMENT: {11/27/2006 7:22:37 PM}	this->AssertValid();
// COMMENT: {11/27/2006 7:22:37 PM}#endif
	return 0;
}

int CSelectedOutput::PushBack(const char* key, const CVar& var)
{
	try
	{
		// check if key is new
		std::map< std::string, size_t >::iterator find;
		find = this->m_mapHeadingToCol.find(std::string(key));
		if (find == this->m_mapHeadingToCol.end())
		{
			// new key(column)
			//
			this->m_mapHeadingToCol.insert(std::map< std::string, size_t >::value_type(std::string(key),
				this->m_mapHeadingToCol.size()) );

			// add heading
			//
			this->m_vecVarHeadings.push_back(CVar(key));


			// add new vector(col)
			//
			this->m_arrayVar.resize(this->m_arrayVar.size() + 1);
			this->m_arrayVar.back().reserve(RESERVE_ROWS);

			// add empty rows if nec
			if (this->m_nRowCount)
				this->m_arrayVar.back().resize(this->m_nRowCount);

			this->m_arrayVar.back().push_back(var);
		}
		else
		{
			if (this->m_arrayVar[find->second].size() == this->m_nRowCount) {
				this->m_arrayVar[find->second].push_back(var);
			}
			else {
				ASSERT(this->m_arrayVar[find->second].size() == this->m_nRowCount + 1);
				this->m_arrayVar[find->second].at(this->m_nRowCount) = var;
			}
		}
		return 0;
	}
	catch(...)
	{
		ASSERT(false);
	}
	return 1; // error
}


int CSelectedOutput::PushBackDouble(const char* key, double value)
{
	CVar v(value);
	return this->PushBack(key, v);
}

int CSelectedOutput::PushBackLong(const char* key, long value)
{
	CVar v(value);
	return this->PushBack(key, v);
}

int CSelectedOutput::PushBackString(const char* key, const char* value)
{
	CVar v(value);
	return this->PushBack(key, v);
}

int CSelectedOutput::PushBackEmpty(const char* key)
{
	CVar v;
	return this->PushBack(key, v);
}

#if defined(_DEBUG)
void CSelectedOutput::Dump(const char* heading)
{
	::OutputDebugStringA(heading);
	std::ostringstream oss;
	oss << *this;
	std::istringstream iss(oss.str());
	std::string line;
	while (std::getline(iss, line)) {
		::OutputDebugStringA(line.c_str());
		::OutputDebugStringA("\n");
	}
}

void CSelectedOutput::AssertValid(void)const
{
	if (size_t cols = this->GetColCount())
	{
		size_t rows = this->m_arrayVar[0].size();
		for (size_t col = 0; col < cols; ++col)
		{
			ASSERT(rows == this->m_arrayVar[col].size());
		}
	}
}
#endif

std::ostream& operator<< (std::ostream &os, const CSelectedOutput &a)
{
#if defined(_WIN32)
	os << "CSelectedOutput(rows=" << a.GetRowCount() << ", cols=" << a.GetColCount() << ")\n";
#endif

	CVar v;
	for (size_t r = 0; r < a.GetRowCount(); ++r) {
		for (size_t c = 0; c < a.GetColCount(); ++c) {
			a.Get((int)r, (int)c, &v);
			os << v << ", ";
			::VarClear(&v);
		}
		os << "\n";
	}
	os << "\n";
	return os;
}
