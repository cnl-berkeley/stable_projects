For detailed instructions on how to use the code provided and how the data is jittered to be in the format required of the analyses carried out in this study, refer closely to comments in each of the appended scripts. If you have any questions or need clarification, please contact William Hastings (will.law.hast@berkeley.edu).

The variable names for each of the files used in the analysis are as follows:

3Factor_impulsivity.csv : contains behavioral measures of ERI. For the purpose of this study, only, 'ID', 'Factor_FeelingsTrigger', 'Factor_LackFollowThrough', and 'Factor_PervInf' are used, which were previously calculated (values range 1-5) based on the responses to questions coded in columns 2:106.
	ID: participant identification
	Factor_FeelingsTrigger: Feelings Trigger Actions
	Factor_LackFollowThrough: Lack of Follow Through
	Factor_PervInf: Pervasive Influence of Feelings
	
urgency_type.csv : contains the number of each variable sulcus in each hemisphere as well as the sulcogyral type.
	sub: participant identification
	hemi: hemisphere
		lh: left hemisphere
		rh: right hemisphere
	type: sulcogyral type
	ios: intermediate orbital sulcus
	pos: posterior orbital sulcus
	sf: sulcus fragmentosis
	nak: the total number of variable sulci excluding sf components
	tot: the total number of variable sulci including sf components

urgency.types.csv : contains the percent occurrence of each sulcogyral type across participants in both, left, or right hemispheres for the purpose of comparison to Nakamura's meta-analysis
	Hemi: Hemisphere
	Type: Sulcogyral type
	Percent: percent occurrence in of each type in corresponding hemisphere(s)

urgency_metrics_new.csv : contains the sulcal metrics of each sulcus. For the purpose of this study, only sub, hemi, cortical_thickness_mean, and sulcal_depth_mean_pct were used.
	sub: participant identification
	hemi: hemisphere
	cortical_thickness_mean: the mean cortical thickness of a sulcus in millimeters
	sulcal_depth_mean_pct: a standardized measure of mean depth of each sulcus in FreeSurfer Units
	olfs: olfactory sulcus
	tolfs: transverse olfactory sulcus
	mosa: anterior portion of the medial orbital sulcus
	mosp: posterior portion of the medial orbital sulcus
	losa: anterior portion of the lateral orbital sulcus
	losp: posterior portion of the lateral orbital sulcus
	tos: transverse orbital sulcus
	ios: intermediate orbital sulcus
	pos: posterior orbital sulcus
	sf: sulcus fragmentosis

urgency_metrics_colfs : identical to urgency_metrics_new.csv, but replaces tolfs and olfs with colfs.
	colfs: olfactory sulcal complex

scid_approach.csv : a csv containing diagnostic and demographic data for each participant. Included also is a codebook used to understand the values listed.
	sub: participant identification
	Education: years participant spent being educated
	Cur_Meds: whether or not the participant currently takes meds
	Lifetime_psychosis: if the participant has a diagnosis of Lifetime Psychosis (not used)
	lifetimeSUD: if the participant has a diagnosis of lifetime substance use disorder
	lifetimeAUD: if the participant has a diagnosis of lifetime alcohol use disorder
	lifetime_MDE_dx: if the participant has a diagnosis of lifetime major depressive episodes	
	Lifetime_anxiety_dx: if the participant has a diagnosis of lifetime anxiety	
	lifetime_sustance_dx: if the participant has a diagnosis of lifetime substance use diagnosis
	lifetime_OC_hoarding_dx: if the participant has a diagnosis of hoarding
	SHEEHAN_MAX: what participant is rated on the Sheehan disability scale
	Gender: the gender the participant most closely identifies with
	Gender_2: string data allowing participant to elaborate on gender identity
	Hispanic: whether the participant identifies as Hispanic
	Race: what race the participant most closely identifies with
	Race_5_TEXT: string data allowing participant to elaborate on race
	Age: participant age in years
	
SCID_codebook.xlsx : a spreadsheet containing the labels associated with each variable value for each variable listed in scid_aproach.csv

MedData_approach.csv : a csv containing information pertaining to what medications the participant is taking or has taken in the past. For the purpose of this study, only a binary "Yes" or "No" was considered, not start date or dosage.
	sub: participant identification
	SSRI_YN: whether or not the participant takes selective serotonin reuptake inhibitors
	SNRI_YN: whether or not the participant takes serotonin-norepinephrine reuptake inhibitors
	Tricyclic_YN: whether or not the participant takes tricyclic antidepressants
	Other_antidep_YN: whether or not the participant takes other antidepressants
	benzodiaz_YN: whether or not the participant takes benzodiazapines
	atypical_antipsych_YN: whether or not the participant takes atypical antipsychotics
	antipsych_YN: whether or not the participant takes antipsychotics
	antiseisureantiepilepticmeds_YN: whether or not the participant takes anti-seizure or anti-epileptic medications
	blood_pressure_YN: whether or not the participant takes blood pressure medication
	statins_YN: whether or not the participant takes statins
	stimulants_YN: whether or not the participant takes stimulants
	oral_contrapceptives_YN: whether or not the participant takes oral contraceptives