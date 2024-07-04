#include "PrecompiledHeaders.h"
#include "WeftKnitWizard.h"
#include "WindowIDs.h"

#include "LoopGeometry.xpm"

using namespace TexGen;

BEGIN_EVENT_TABLE(CWeftKnitWizard, wxWizard)
EVT_WIZARD_PAGE_CHANGING(wxID_ANY, CWeftKnitWizard::OnWizardPageChanging)
EVT_INIT_DIALOG(CWeftKnitWizard::OnInit)
EVT_CHECKBOX(ID_Refine, CWeftKnitWizard::OnRefine)
EVT_CHECKBOX(ID_DefaultDomain, CWeftKnitWizard::OnDomain)
EVT_TEXT(ID_WaleHeight, CWeftKnitWizard::OnWaleHeightChanged)
EVT_TEXT(ID_CourseWidth, CWeftKnitWizard::OnCourseWidthChanged)
EVT_TEXT(ID_LoopHeight, CWeftKnitWizard::OnLoopHeightChanged)
EVT_TEXT(ID_Thickness, CWeftKnitWizard::OnThicknessChanged)
END_EVENT_TABLE()

CWeftKnitWizard::CWeftKnitWizard(wxWindow* parent, wxWindowID id)
	: wxWizard(parent, id, wxT("Weft Knit Wizard"), wxBitmap(LoopGeometry_xpm))
	, m_bCreateDomain(true)
	, m_pFirstPage(NULL)
	, m_WaleHeight(wxT("1"))
	, m_LoopHeight(wxT("1.2"))
	, m_CourseWidth(wxT("1"))
	, m_YarnThickness(wxT("0.2"))
	, m_GapSize(wxT("0"))
	, m_LoopModel(RAVANDI_2021)
	, m_bRefine(true)
	, m_bWaleHeightChanged(false)
	, m_bCourseWidthChanged(false)
	, m_bLoopHeightChanged(false)
	, m_bThicknessChanged(false)
	, m_bBuiltPages(false)
{
	BuildPages();
	GetPageAreaSizer()->Add(m_pFirstPage);
}

CWeftKnitWizard::~CWeftKnitWizard(void)
{
	
}

void CWeftKnitWizard::OnWizardPageChanging(wxWizardEvent& event)
{
	
}

bool CWeftKnitWizard::HasNextPage(wxWizardPage *page)
{
	if (page == m_pFirstPage)
		return true;
	return wxWizard::HasNextPage(page);
}


void CWeftKnitWizard::BuildPages()
{
	m_pFirstPage = BuildFirstPage();
	m_pSecondPage = BuildSecondPage();

	wxWizardPageSimple::Chain(m_pFirstPage, m_pSecondPage);

	m_bBuiltPages = true;
}


wxWizardPageSimple* CWeftKnitWizard::BuildFirstPage()
{
	wxWizardPageSimple *pPage = new wxWizardPageSimple(this);

	wxBoxSizer *pMainSizer = new wxBoxSizer(wxVERTICAL);
	wxSizerFlags SizerFlags(0);

	SizerFlags.Border();
	SizerFlags.Expand();

	wxSizer *pSubSizer;

	pMainSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Specify loop model and loop quantities for the weft knit textile.")), SizerFlags);

	SizerFlags.Align(wxALIGN_CENTER_VERTICAL);

	{
		wxArrayString loopModelChoices;
		loopModelChoices.Add(wxT("Ravandi et al. 2022"));
		
		pMainSizer->Add(m_pLoopModelRadio = new wxRadioBox(pPage, ID_LoopModel, wxT("Loop Model"), wxDefaultPosition, wxDefaultSize, loopModelChoices, 2, wxRA_SPECIFY_COLS), SizerFlags);
		m_pLoopModelRadio->SetToolTip(wxT("Specifies which loop model is used to determine master node points."));

	}
	
	{
		pSubSizer = new wxFlexGridSizer(2);

		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Wales:")), SizerFlags);
		pSubSizer->Add(m_pWalesSpin = new wxSpinCtrl(pPage, wxID_ANY, wxT(""), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 100), SizerFlags);
		m_pWalesSpin->SetToolTip(wxT("Controls the number of wales contained in the unit cell."));

		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Courses:")), SizerFlags);
		pSubSizer->Add(m_pCoursesSpin = new wxSpinCtrl(pPage, wxID_ANY, wxT(""), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 100), SizerFlags);
		m_pCoursesSpin->SetToolTip(wxT("Controls the number of courses contained in the unit cell."));

		wxTextCtrl* pControl;

		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Wale Height:")), SizerFlags);
		pSubSizer->Add(pControl = new wxTextCtrl(pPage, ID_Spacing, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC, &m_WaleHeight)), SizerFlags);
		pControl->SetToolTip(wxT("Sets the height of one wale."));

		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Loop Height:")), SizerFlags);
		pSubSizer->Add(pControl = new wxTextCtrl(pPage, ID_Width, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC, &m_LoopHeight)), SizerFlags);
		pControl->SetToolTip(wxT("Sets the height of one loop."));

		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Course Width:")), SizerFlags);
		pSubSizer->Add(pControl = new wxTextCtrl(pPage, ID_Width, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC, &m_CourseWidth)), SizerFlags);
		pControl->SetToolTip(wxT("Sets the width of one course."));

		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Yarn Thickness:")), SizerFlags);
		pSubSizer->Add(pControl = new wxTextCtrl(pPage, ID_Thickness, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC, &m_YarnThickness)), SizerFlags);
		pControl->SetToolTip(wxT("Sets the thickness of the yarn."));

		pMainSizer->Add(pSubSizer, SizerFlags);
	}

	SizerFlags.Align(0);

	SizerFlags.Align(wxALIGN_CENTER_VERTICAL);
	
	{
		pSubSizer = new wxFlexGridSizer(3);

		wxCheckBox* pDomainBox;
		pSubSizer->Add(pDomainBox = new wxCheckBox(pPage, ID_DefaultDomain, wxT("Create default domain"), wxDefaultPosition, wxDefaultSize, 0, wxGenericValidator(&m_bCreateDomain)), SizerFlags);
		pDomainBox->SetToolTip(wxT("Assign the default domain to the textile. The default domain is the same size as the unit cell in the x and y directions and and corresponds to the fabric thickness in the z direction."));

		pSubSizer->AddSpacer(0);
		pSubSizer->AddSpacer(0);

		/*wxCheckBox* pRefineBox;
		pSubSizer->Add(pRefineBox = new wxCheckBox(pPage, ID_Refine, wxT("Refine model"), wxDefaultPosition, wxDefaultSize, 0, wxGenericValidator(&m_bRefine)), SizerFlags);
		pRefineBox->SetToolTip(wxT("Refine the model cross sections such that yarn volumes don't intersect.\nThe resulting geometry always contains a minimum gap size between two yarns."));
		
		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Gap size:")), SizerFlags);
		wxTextCtrl* pGapSize;
		pSubSizer->Add(pGapSize = new wxTextCtrl(pPage, ID_GapSize, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC, &m_GapSize)), SizerFlags);
		pGapSize->SetToolTip(wxT("Adjusts the minimum gap size between two yarns when refine model is enabled."));

		if (!m_bRefine)
			pGapSize->Disable();*/
	}

	pMainSizer->Add(pSubSizer, SizerFlags);
	SizerFlags.Align(0);

	pPage->SetSizer(pMainSizer);
	pMainSizer->Fit(pPage);

	m_pWalesSpin->SetValue(2);
	m_pCoursesSpin->SetValue(2);

	return pPage;
}

wxWizardPageSimple* CWeftKnitWizard::BuildSecondPage()
{
	wxWizardPageSimple *pPage = new wxWizardPageSimple(this);

	wxBoxSizer *pMainSizer = new wxBoxSizer(wxVERTICAL);
	wxSizerFlags SizerFlags(0);

	SizerFlags.Border();
	SizerFlags.Expand();

	wxSizer *pSubSizer;

	pMainSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Specify rendered mesh resolution.")), SizerFlags);
	SizerFlags.Align(wxALIGN_CENTER_VERTICAL);
	{
		pSubSizer = new wxFlexGridSizer(2);

		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Number of slave nodes:")), SizerFlags);
		pSubSizer->Add(m_pNumSlaveNodesSpin = new wxSpinCtrl(pPage, wxID_ANY, wxT(""), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 100), SizerFlags);
		m_pWalesSpin->SetToolTip(wxT("Controls the number of slave nodes to add in between master nodes."));

		pSubSizer->Add(new wxStaticText(pPage, wxID_ANY, wxT("Number of section points:")), SizerFlags);
		pSubSizer->Add(m_pNumSectionPointsSpin = new wxSpinCtrl(pPage, wxID_ANY, wxT(""), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 100), SizerFlags);
		m_pCoursesSpin->SetToolTip(wxT("Controls the number of section points to use across the cross-section of the yarn."));

		pMainSizer->Add(pSubSizer, SizerFlags);
	}

	pPage->SetSizer(pMainSizer);
	pMainSizer->Fit(pPage);

	m_pNumSlaveNodesSpin->SetValue(50);
	m_pNumSectionPointsSpin->SetValue(20);

	return pPage;
}

bool CWeftKnitWizard::RunIt()
{
	return RunWizard(m_pFirstPage);
}

string CWeftKnitWizard::GetCreateTextileCommand(string ExistingTextile)
{
	stringstream StringStream;

	int iWales = m_pWalesSpin->GetValue();
	int iCourses = m_pCoursesSpin->GetValue();
	double dWaleHeight, dLoopHeight, dCourseWidth, dThickness;
	m_WaleHeight.ToDouble(&dWaleHeight);
	m_LoopHeight.ToDouble(&dLoopHeight);
	m_CourseWidth.ToDouble(&dCourseWidth);
	m_YarnThickness.ToDouble(&dThickness);
	int iNumSlaveNodes = m_pNumSlaveNodesSpin->GetValue();
	int iNumSectionPoints = m_pNumSectionPointsSpin->GetValue();
	
	StringStream << "WeftKnit = CTextileWeftKnit(" << iWales << "," << iCourses << "," << dWaleHeight << "," << dLoopHeight << "," << dCourseWidth << "," << dThickness << ")" << endl;

	StringStream << "WeftKnit.SetNumSlaveNodes(" << iNumSlaveNodes << ")" << endl;

	StringStream << "WeftKnit.SetNumSectionPoints(" << iNumSectionPoints << ")" << endl;

	if (m_bCreateDomain)
	{
		StringStream << "WeftKnit.AssignDefaultDomain()" << endl;
	}

	/*if (m_bRefine)
	{
		StringStream << "WeftKnit.RefineTextile()" << endl;
	}*/

	switch (m_pLoopModelRadio->GetSelection())
	{
		case RAVANDI_2021:
		{
			StringStream << "WeftKnit.SetLoopModel(RAVANDI_2021)" << endl;
		}
	}
	

	StringStream << "AddTextile(WeftKnit)" << endl;

	return StringStream.str();
}

void CWeftKnitWizard::LoadSettings(const CTextileWeftKnit &WeftKnit)
{

}



void CWeftKnitWizard::OnRefine(wxCommandEvent& event)
{
	RefreshGapTextBox();
}


void CWeftKnitWizard::OnDomain(wxCommandEvent& event)
{

}

void CWeftKnitWizard::RefreshGapTextBox()
{
	wxCheckBox* pRefine = (wxCheckBox*)FindWindowById(ID_Refine, this);
	wxTextCtrl* pGapSize = (wxTextCtrl*)FindWindowById(ID_GapSize, this);

	if (pRefine)
	{
		pRefine->Enable();

		if (pRefine->GetValue())
			pGapSize->Enable();
		else
			pGapSize->Disable();
	}
}