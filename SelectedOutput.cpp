// SelectedOutput.cpp: implementation of the CSelectedOutput class.
//
//////////////////////////////////////////////////////////////////////
#if defined(_DEBUG)
#pragma warning(disable : 4786) // disable truncation warning
#include <sstream>              // std::ostringstream
#include <windows.h>            // ::OutputDebugString
#endif

#include <stdarg.h>
#include <stdio.h>

#include "SelectedOutput.hxx"
#include "phreeqcns.hxx"

const size_t RESERVE_ROWS = 80;
const size_t RESERVE_COLS = 80;


int EndRow(void);
void AddSelectedOutput(const char* name, const char* format, va_list argptr);
int warning_msg (const char *err_str);

int EndRow(void)
{
	if (CSelectedOutput::Instance()->GetRowCount() <= 1) {
		// ensure all user_punch headings are included
		for (int i = n_user_punch_index; i < user_punch_count_headings; ++i) {
			CSelectedOutput::Instance()->PushBackEmpty(user_punch_headings[i]);
		}
	}
	return CSelectedOutput::Instance()->EndRow();
}

int PushBackDouble(const char* key, double dVal)
{
	return CSelectedOutput::Instance()->PushBackDouble(key, dVal);
}

int PushBackLong(const char* key, long lVal)
{
	return CSelectedOutput::Instance()->PushBackLong(key, lVal);
}

int PushBackString(const char* key, const char* sVal)
{
	return CSelectedOutput::Instance()->PushBackString(key, sVal);
}

int PushBackEmpty(const char* key)
{
	return CSelectedOutput::Instance()->PushBackEmpty(key);
}


void AddSelectedOutput(const char* name, const char* format, va_list argptr)
{
	int bInt;
	int bDouble;
	int bString;

	int state;
	int bLongDouble;
	char ch;


	/* state values
	0 Haven't found start(%)
	1 Just read start(%)
	2 Just read Flags(-0+ #) (zero or more)
	3 Just read Width
	4 Just read Precision start (.)
	5 Just read Size modifier
	6 Just read Type
	*/

	if (name == NULL) {
		return;
	}

	bDouble = 0;
	bInt = 0;
	bString = 0;

	bLongDouble = 0;

	state = 0;
	ch = *format++;
	while (ch != '\0') {
		switch (state) {
	case 0: /* looking for Start specification (%) */
		switch (ch) {
	case '%':
		state = 1;
		break;
	default:
		break;
		}
		ch = *format++;
		break;
	case 1: /* reading Flags (zero or more(-,+,0,# or space)) */
		switch (ch) {
	case '-': case '0': case '+': case ' ': case '#':
		ch = *format++;
		break;
	default:
		state = 2;
		break;
		}
		break;
	case 2: /* reading Minimum field width (decimal integer constant) */
		switch (ch) {
	case '.':
		state = 3;
		ch = *format++;
		break;
	case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
		ch = *format++;
		break;
	default:
		state = 4;
		break;
		}
		break;
	case 3: /* reading Precision specification (period already read) */
		switch (ch)
		{
		case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
			ch = *format++;
			break;
		default:
			state = 4;
			break;
		}
		break;
	case 4: /* reading Size modifier */
		switch (ch)
		{
		case 'l':
			ch = *format++;
			break;
		case 'L':
			bLongDouble = 1;
			ch = *format++;
			break;
		case 'h':
			ch = *format++;
			break;
		}
		state = 5;
		break;
	case 5: /* reading Conversion letter */
		switch (ch) {
	case 'c':
		break;
	case 'd':
	case 'i':
		bInt = 1;
		break;
	case 'n':
	case 'o':
	case 'p':
		break;
	case 's':
		bString = 1;
		break;
	case 'u':
	case 'x':
	case 'X':
	case '%':
		break;
	case 'f':
	case 'e':
	case 'E':
	case 'g':
	case 'G':
		bDouble = 1;
		break;
	default:
		ASSERT(false);
		break;
		}
		ch = '\0';  /* done */
		break;
		}
	}

	if (bDouble) {
		double valDouble;

		if (bLongDouble) {
			valDouble = (double)va_arg(argptr, long double);
		}
		else {
			valDouble = va_arg(argptr, double);
		}

		CSelectedOutput::Instance()->PushBackDouble(name, valDouble);
	}
	else if (bInt) {
		int valInt;
		valInt = va_arg(argptr, int);

		CSelectedOutput::Instance()->PushBackLong(name, (long)valInt);
	}
	else if (bString) {
		char* valString;
		valString = (char *)va_arg(argptr, char *);

		CSelectedOutput::Instance()->PushBackString(name, valString);
	}
	else {
		ASSERT(false);
		CSelectedOutput::Instance()->PushBackEmpty(name);
	}
}

// COMMENT: {11/16/2004 10:18:22 PM}CSelectedOutput CSelectedOutput::singleton;
CSelectedOutput* CSelectedOutput::s_instance = 0;

CSelectedOutput* CSelectedOutput::Instance()
{
	if (s_instance == 0)
	{
		s_instance = new CSelectedOutput;
	}
	return s_instance;
}

void CSelectedOutput::Release()
{
	if (s_instance)
	{
		delete s_instance;
		s_instance = 0;
	}
}

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
