package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.APPROVEDCOUNTRY;
import static com.hartwig.hmftools.vicc.database.Tables.ASSOCIATION;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAPART1;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAPART2;
import static com.hartwig.hmftools.vicc.database.Tables.CGI;
import static com.hartwig.hmftools.vicc.database.Tables.CGICDNA;
import static com.hartwig.hmftools.vicc.database.Tables.CGIGDNA;
import static com.hartwig.hmftools.vicc.database.Tables.CGIINDIVIDUALMUTATION;
import static com.hartwig.hmftools.vicc.database.Tables.CGIINFO;
import static com.hartwig.hmftools.vicc.database.Tables.CGIREGION;
import static com.hartwig.hmftools.vicc.database.Tables.CGISTRAND;
import static com.hartwig.hmftools.vicc.database.Tables.CGITRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.CIVIC;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICASSERTIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCLINICALTRIAL;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCLINVARENTRIES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDESCRIPTION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDISEASE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDRUGS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICERROR;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSCLINICALTRIAL;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSPUBLICATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICHGVSEXPRESSIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDON;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIED;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWED;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLIFECYCLEACTIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICPUBLICATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTALIASES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSTYPES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSVARIANTS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTTYPES;
import static com.hartwig.hmftools.vicc.database.Tables.DEVTAG;
import static com.hartwig.hmftools.vicc.database.Tables.ENVIRONMENTALCONTEXT;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCE;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCEINFO;
import static com.hartwig.hmftools.vicc.database.Tables.EVIDENCETYPE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURE;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURENAME;
import static com.hartwig.hmftools.vicc.database.Tables.GENE;
import static com.hartwig.hmftools.vicc.database.Tables.GENEIDENTIFIER;
import static com.hartwig.hmftools.vicc.database.Tables.HIERARCHY;
import static com.hartwig.hmftools.vicc.database.Tables.JAX;
import static com.hartwig.hmftools.vicc.database.Tables.JAXINDICATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXREFERENCES;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTHERAPY;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSINDICATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSTHERAPIES;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSVARIANTREQUIREMENTDETAILS;
import static com.hartwig.hmftools.vicc.database.Tables.LINK;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCH;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONALT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONCHR;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONCOSMICID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONDBSNP;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEXON;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEXONICFUNC;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPARENTS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPATHOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPOPFREQMAX;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONREF;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONSOURCES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONSTART;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCRITERIAMET;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCRITERIAUNMET;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHEXTERNALID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDECONDITION0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDECONDITION1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEDRUG1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEGENE0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEMUTATION0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEMUTATION1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDESTAGE0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINSTUTITION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSCDNA;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONSINFO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONABAND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONACHR;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONACOORD;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONAGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONAORI;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONAREG;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONATX;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBBAND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBCHR;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBCOORD;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBORI;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBREG;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBTX;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONINS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSGRCH37LOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSMUTATIONTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSPARENTS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSPATHOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSSYNONYMS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATA;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATAFULLAA;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATAGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSMAP;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHPREFELANCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTAGS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTHERAPEUTICCONTEXT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTIEREXPLANATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSALTERATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSCONTACT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSGEO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSINTERVATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSLOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSLOCATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOTHERGROUPLABEL;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOTHERNAME;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOVERALLCONTACT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSTAGS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOCONSEQUENCES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOFUSIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOLOCATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKB;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCONSEQUENCESBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCONSEQUENCESCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBDRUGABSTRACTSCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEALIASESBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEALIASESCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENECLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBVARIANTBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBVARIANTCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PHENOTYPETYPE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKB;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBGENE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBTISSUE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBTUMOR;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBVARIANT;
import static com.hartwig.hmftools.vicc.database.Tables.PROVENANCE;
import static com.hartwig.hmftools.vicc.database.Tables.PUBLICATIONURL;
import static com.hartwig.hmftools.vicc.database.Tables.SAGE;
import static com.hartwig.hmftools.vicc.database.Tables.SEQUENCEONTOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.SYNONYM;
import static com.hartwig.hmftools.vicc.database.Tables.TAG;
import static com.hartwig.hmftools.vicc.database.Tables.TAXONOMY;
import static com.hartwig.hmftools.vicc.database.Tables.VICCENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.vicc.ViccJsonToSQLImporter;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.BRCA;
import com.hartwig.hmftools.vicc.datamodel.Cgi;
import com.hartwig.hmftools.vicc.datamodel.Civic;
import com.hartwig.hmftools.vicc.datamodel.CivicClinicalTrial;
import com.hartwig.hmftools.vicc.datamodel.CivicDrugs;
import com.hartwig.hmftools.vicc.datamodel.CivicEvidenceItems;
import com.hartwig.hmftools.vicc.datamodel.CivicSource;
import com.hartwig.hmftools.vicc.datamodel.CivicUser;
import com.hartwig.hmftools.vicc.datamodel.CivicVariantGroup;
import com.hartwig.hmftools.vicc.datamodel.CivicVariantTypes;
import com.hartwig.hmftools.vicc.datamodel.CivicVariants;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.Jax;
import com.hartwig.hmftools.vicc.datamodel.JaxReferences;
import com.hartwig.hmftools.vicc.datamodel.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsIndications;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsTherapies;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsVariantRequirementDetails;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAreg;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchBreg;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchExonInfo;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchFusions;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchGRch37Location;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchLocations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchMutations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchParents;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTags;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTranscriptConsequencesGRCH37;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsIntervation;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsLocations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsTags;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchWGSaMap;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchWGSadataLocation;
import com.hartwig.hmftools.vicc.datamodel.OncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.OncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.OncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.OncoKbDrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.Oncokb;
import com.hartwig.hmftools.vicc.datamodel.OncokbGene;
import com.hartwig.hmftools.vicc.datamodel.OncokbVariant;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.PmkbGene;
import com.hartwig.hmftools.vicc.datamodel.PmkbTissue;
import com.hartwig.hmftools.vicc.datamodel.PmkbTumor;
import com.hartwig.hmftools.vicc.datamodel.PmkbVariant;
import com.hartwig.hmftools.vicc.datamodel.Sage;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;

public class ViccDAO {

    private static final Logger LOGGER = LogManager.getLogger(ViccJsonToSQLImporter.class);

    @NotNull
    private final DSLContext context;

    public static ViccDAO connectToViccDAO(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
            throws SQLException {
        final Connection conn = DriverManager.getConnection(url, userName, password);
        final String catalog = conn.getCatalog();

        LOGGER.debug("Connecting to database {}", catalog);
        return new ViccDAO(DSL.using(conn, SQLDialect.MYSQL));
    }

    private ViccDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void writeViccEntry(@NotNull ViccEntry viccEntry) {
        int id = context.insertInto(VICCENTRY, VICCENTRY.SOURCE)
                .values(viccEntry.source())
                .returning(VICCENTRY.ID)
                .fetchOne()
                .getValue(VICCENTRY.ID);
        LOGGER.info("vicc entry id: " + id);
        writeTags(id, viccEntry.tags());
        writeDevTags(id, viccEntry.devTags());
        writeGeneIdentifiers(id, viccEntry.geneIdentifiers());
        writeGenes(id, viccEntry.genes());
        writeFeatureNames(id, viccEntry.featureNames());
        writeFeatures(id, viccEntry.features());
        writeAssociation(id, viccEntry.association());
        writeKbSpecificObject(id, viccEntry.KbSpecificObject());
    }

    private void importSageinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Sage sage = (Sage) object;
        context.insertInto(SAGE,
                SAGE.ENTREZID,
                SAGE.CLINICALMANIFESTATION,
                SAGE.PUBLICATIONURL,
                SAGE.GERMLINEORSOMATIC,
                SAGE.EVIDENCELABEL,
                SAGE.DRUGLABEL,
                SAGE.RESPONSETYPE,
                SAGE.GENE,
                SAGE.VICCENTRYID)
                .values(sage.entrezId(),
                        sage.clinicalManifestation(),
                        sage.publicationUrl(),
                        sage.germlineOrSomatic(),
                        sage.evidenceLabel(),
                        sage.drugLabels(),
                        sage.responseType(),
                        sage.gene(),
                        viccEntryId)
                .execute();
    }

    private void importBRCAinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        BRCA brca = (BRCA) object;
        context.insertInto(BRCAPART1,
                BRCAPART1.VARIANT_FREQUENCY_LOVD,
                BRCAPART1.ALLELE_FREQUENCY_FIN_EXAC,
                BRCAPART1.CLINVARACCESSION_ENIGMA,
                BRCAPART1.HOMOZYGOUS_COUNT_AFR_EXAC,
                BRCAPART1.BX_ID_EXAC,
                BRCAPART1.VARIANT_IN_LOVD,
                BRCAPART1.ALLELE_FREQUENCY_AFR_EXAC,
                BRCAPART1.DBID_LOVD,
                BRCAPART1.CHR,
                BRCAPART1.BX_ID_ENIGMA,
                BRCAPART1.CO_OCCURRENCE_LR_EXLOVD,
                BRCAPART1.HOMOZYGOUS_COUNT_EAS_EXAC,
                BRCAPART1.SUBMITTER_CLINVAR,
                BRCAPART1.ALLELE_FREQUENCY_EAS_EXAC,
                BRCAPART1.HG37_END,
                BRCAPART1.SUBMITTERS_LOVD,
                BRCAPART1.CLINICAL_CLASSIFICATION_BIC,
                BRCAPART1.HOMOZYGOUS_COUNT_NFE_EXAC,
                BRCAPART1.ALLELE_COUNT_SAS_EXAC,
                BRCAPART1.METHOD_CLINVAR,
                BRCAPART1.ALLELE_COUNT_NFE_EXAC,
                BRCAPART1.PATHOGENICITY_ALL,
                BRCAPART1.GERMLINE_OR_SOMATIC_BIC,
                BRCAPART1.HOMOZYGOUS_COUNT_SAS_EXAC,
                BRCAPART1.BIC_NOMENCLATURE,
                BRCAPART1.ASSERTION_METHOD_ENIGMA,
                BRCAPART1.LITERATURE_SOURCE_EXLOVD,
                BRCAPART1.CHANGE_TYPE_ID,
                BRCAPART1.COLLECTION_METHOD_ENIGMA,
                BRCAPART1.SUM_FAMILY_LR_EXLOVD,
                BRCAPART1.HGVS_CDNA_LOVD,
                BRCAPART1.HOMOZYGOUS_COUNT_FIN_EXAC,
                BRCAPART1.EAS_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART1.ETHNICITY_BIC,
                BRCAPART1.INDIVIDUALS_LOVD,
                BRCAPART1.VARIANT_IN_EXAC,
                BRCAPART1.URL_ENIGMA,
                BRCAPART1.ALLELE_ORIGIN_CLINVAR,
                BRCAPART1.ALLELE_FREQUENCY_AMR_EXAC,
                BRCAPART1.VARIANT_IN_1000_GENOMES,
                BRCAPART1.AFR_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART1.BX_ID_EXLOVD,
                BRCAPART1.SOURCE,
                BRCAPART1.CONDITION_ID_VALUE_ENIGMA,
                BRCAPART1.HGVS_PROTEIN,
                BRCAPART1.REF,
                BRCAPART1.ALLELE_NUMBER_AFR_EXAC,
                BRCAPART1.ALLELE_COUNT_AFR_EXAC,
                BRCAPART1.BX_ID_LOVD,
                BRCAPART1.SYNONYMS,
                BRCAPART1.GENE_SYMBOL,
                BRCAPART1.COMMENT_ON_CLINICAL_SIGNIFICANCE_ENIGMA,
                BRCAPART1.MISSENSE_ANALYSIS_PRIOR_PROBABILITY_EXLOVD,
                BRCAPART1.ALLELE_NUMBER_FIN_EXAC,
                BRCAPART1.POSTERIOR_PROBABILITY_EXLOVD,
                BRCAPART1.POLYPHEN_SCORE,
                BRCAPART1.REFERENCE_SEQUENCE,
                BRCAPART1.ALLELE_COUNT_EAS_EXAC,
                BRCAPART1.HG38_END,
                BRCAPART1.HGVS_CDNA,
                BRCAPART1.FUNCTIONAL_ANALYSIS_TECHNIQUE_LOVD,
                BRCAPART1.SAS_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART1.RNA_LOVD,
                BRCAPART1.COMBINED_PRIOR_PROBABLILITY_EXLOVD,
                BRCAPART1.BX_ID_CLINVAR,
                BRCAPART1.IARC_CLASS_EXLOVD,
                BRCAPART1.BX_ID_BIC,
                BRCAPART1.VICCENTRYID)
                .values(brca.brcApart1().Variant_frequency_LOVD(),
                        brca.brcApart1().Allele_frequency_FIN_ExAC(),
                        brca.brcApart1().ClinVarAccession_ENIGMA(),
                        brca.brcApart1().Homozygous_count_AFR_ExAC(),
                        brca.brcApart1().BX_ID_ExAC(),
                        brca.brcApart1().Variant_in_LOVD(),
                        brca.brcApart1().Allele_frequency_AFR_ExAC(),
                        brca.brcApart2().DBID_LOVD(),
                        brca.brcApart1().Chr(),
                        brca.brcApart1().BX_ID_ENIGMA(),
                        brca.brcApart1().Co_occurrence_LR_exLOVD(),
                        brca.brcApart1().Homozygous_count_EAS_ExAC(),
                        brca.brcApart1().Submitter_ClinVar(),
                        brca.brcApart1().Allele_frequency_EAS_ExAC(),
                        brca.brcApart1().Hg37_End(),
                        brca.brcApart1().Submitters_LOVD(),
                        brca.brcApart1().Clinical_classification_BIC(),
                        brca.brcApart1().Homozygous_count_NFE_ExAC(),
                        brca.brcApart1().Allele_count_SAS_ExAC(),
                        brca.brcApart1().Method_ClinVar(),
                        brca.brcApart1().Allele_count_NFE_ExAC(),
                        brca.brcApart1().Pathogenicity_all(),
                        brca.brcApart1().Germline_or_Somatic_BIC(),
                        brca.brcApart1().Homozygous_count_SAS_ExAC(),
                        brca.brcApart1().BIC_Nomenclature(),
                        brca.brcApart1().Assertion_method_ENIGMA(),
                        brca.brcApart1().Literature_source_exLOVD(),
                        brca.brcApart1().Change_Type_id(),
                        brca.brcApart1().Collection_method_ENIGMA(),
                        brca.brcApart1().Sum_family_LR_exLOVD(),
                        brca.brcApart1().HGVS_cDNA_LOVD(),
                        brca.brcApart1().Homozygous_count_FIN_ExAC(),
                        brca.brcApart1().EAS_Allele_frequency_1000_Genomes(),
                        brca.brcApart1().Ethnicity_BIC(),
                        brca.brcApart1().Individuals_LOVD(),
                        brca.brcApart1().Variant_in_ExAC(),
                        brca.brcApart1().URL_ENIGMA(),
                        brca.brcApart1().Allele_Origin_ClinVar(),
                        brca.brcApart1().Allele_frequency_AMR_ExAC(),
                        brca.brcApart1().Variant_in_1000_Genomes(),
                        brca.brcApart1().AFR_Allele_frequency_1000_Genomes(),
                        brca.brcApart1().BX_ID_exLOVD(),
                        brca.brcApart1().Source(),
                        brca.brcApart1().Condition_ID_value_ENIGMA(),
                        brca.brcApart1().HGVS_Protein(),
                        brca.brcApart1().Ref(),
                        brca.brcApart1().Allele_number_AFR_ExAC(),
                        brca.brcApart1().Allele_count_AFR_ExAC(),
                        brca.brcApart1().BX_ID_LOVD(),
                        brca.brcApart1().Synonyms(),
                        brca.brcApart1().Gene_Symbol(),
                        brca.brcApart1().Comment_on_clinical_significance_ENIGMA(),
                        brca.brcApart1().Missense_analysis_prior_probability_exLOVD(),
                        brca.brcApart1().Allele_number_FIN_ExAC(),
                        brca.brcApart1().Posterior_probability_exLOVD(),
                        brca.brcApart1().Polyphen_Score(),
                        brca.brcApart1().Reference_Sequence(),
                        brca.brcApart1().Allele_count_EAS_ExAC(),
                        brca.brcApart1().Hg38_End(),
                        brca.brcApart1().HGVS_cDNA(),
                        brca.brcApart1().Functional_analysis_technique_LOVD(),
                        brca.brcApart1().SAS_Allele_frequency_1000_Genomes(),
                        brca.brcApart1().RNA_LOVD(),
                        brca.brcApart1().Combined_prior_probablility_exLOVD(),
                        brca.brcApart1().BX_ID_ClinVar(),
                        brca.brcApart1().IARC_class_exLOVD(),
                        brca.brcApart1().BX_ID_BIC(),
                        viccEntryId)
                .execute();

        context.insertInto(BRCAPART2,
                BRCAPART2.SIFT_PREDICTION,
                BRCAPART2.ALLELE_NUMBER_NFE_EXAC,
                BRCAPART2.ALLELE_ORIGIN_ENIGMA,
                BRCAPART2.ALLELE_NUMBER_OTH_EXAC,
                BRCAPART2.HG36_END,
                BRCAPART2.ALLELE_FREQUENCY_SAS_EXAC,
                BRCAPART2.DATE_LAST_UPDATED_CLINVAR,
                BRCAPART2.ALLELE_NUMBER_EAS_EXAC,
                BRCAPART2.ALLELE_FREQUENCY_OTH_EXAC,
                BRCAPART2.SOURCE_URL,
                BRCAPART2.SCV_CLINVAR,
                BRCAPART2.PATHOGENICITY_EXPERT,
                BRCAPART2.ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART2.FUNCTIONAL_ANALYSIS_RESULT_LOVD,
                BRCAPART2.AMR_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART2.VARIANT_IN_ESP,
                BRCAPART2.VARIANT_IN_BIC,
                BRCAPART2.CLINICAL_SIGNIFICANCE_ENIGMA,
                BRCAPART2.MAX_ALLELE_FREQUENCY,
                BRCAPART2.ALLELE_COUNT_AMR_EXAC,
                BRCAPART2.VARIANT_IN_ENIGMA,
                BRCAPART2.BX_ID_ESP,
                BRCAPART2.PATIENT_NATIONALITY_BIC,
                BRCAPART2.BX_ID_1000_GENOMES,
                BRCAPART2.GENOMIC_COORDINATE_HG37,
                BRCAPART2.GENOMIC_COORDINATE_HG36,
                BRCAPART2.EUR_ALLELE_FREQUENCY_1000_GENOMES,
                BRCAPART2.NUMBER_OF_FAMILY_MEMBER_CARRYING_MUTATION_BIC,
                BRCAPART2.SEGREGATION_LR_EXLOVD,
                BRCAPART2.ALLELE_FREQUENCY,
                BRCAPART2.MINOR_ALLELE_FREQUENCY_PERCENT_ESP,
                BRCAPART2.ALLELE_FREQUENCY_EXAC,
                BRCAPART2.MUTATION_TYPE_BIC,
                BRCAPART2.ASSERTION_METHOD_CITATION_ENIGMA,
                BRCAPART2.CONDITION_ID_TYPE_ENIGMA,
                BRCAPART2.ALLELE_COUNT_OTH_EXAC,
                BRCAPART2.HGVS_PROTEIN_LOVD,
                BRCAPART2.VARIANT_IN_CLINVAR,
                BRCAPART2.CLINICAL_IMPORTANCE_BIC,
                BRCAPART2.DISCORDANT,
                BRCAPART2.ALLELE_COUNT_FIN_EXAC,
                BRCAPART2.CONDITION_CATEGORY_ENIGMA,
                BRCAPART2.ALLELE_FREQUENCY_ESP,
                BRCAPART2.HOMOZYGOUS_COUNT_OTH_EXAC,
                BRCAPART2.GENETIC_ORIGIN_LOVD,
                BRCAPART2.HOMOZYGOUS_COUNT_AMR_EXAC,
                BRCAPART2.CLINICAL_SIGNIFICANCE_CLINVAR,
                BRCAPART2.AA_ALLELE_FREQUENCY_ESP,
                BRCAPART2.PROTEIN_CHANGE,
                BRCAPART2.VARIANT_IN_EXLOVD,
                BRCAPART2.EA_ALLELE_FREQUENCY_ESP,
                BRCAPART2.HGVS_RNA,
                BRCAPART2.CLINICAL_SIGNIFICANCE_CITATIONS_ENIGMA,
                BRCAPART2.VARIANT_EFFECT_LOVD,
                BRCAPART2.POLYPHEN_PREDICTION,
                BRCAPART2.DATA_RELEASE_ID,
                BRCAPART2.HG37_START,
                BRCAPART2.HG36_START,
                BRCAPART2.SIFT_SCORE,
                BRCAPART2.GENOMIC_COORDINATE_HG38,
                BRCAPART2.ALT,
                BRCAPART2.LITERATURE_CITATION_BIC,
                BRCAPART2.VARIANT_HAPLOTYPE_LOVD,
                BRCAPART2.ALLELE_FREQUENCY_NFE_EXAC,
                BRCAPART2.HG38_START,
                BRCAPART2.POS,
                BRCAPART2.DATE_LAST_EVALUATED_ENIGMA,
                BRCAPART2.ALLELE_NUMBER_SAS_EXAC,
                BRCAPART2.ALLELE_NUMBER_AMR_EXAC,
                BRCAPART2.VICCENTRYID)
                .values(brca.brcApart1().Sift_Prediction(),
                        brca.brcApart1().Allele_number_NFE_ExAC(),
                        brca.brcApart1().Allele_origin_ENIGMA(),
                        brca.brcApart1().Allele_number_OTH_ExAC(),
                        brca.brcApart1().Hg36_End(),
                        brca.brcApart1().Allele_frequency_SAS_ExAC(),
                        brca.brcApart1().Date_Last_Updated_ClinVar(),
                        brca.brcApart1().Allele_number_EAS_ExAC(),
                        brca.brcApart1().Allele_frequency_OTH_ExAC(),
                        brca.brcApart1().Source_URL(),
                        brca.brcApart1().SCV_ClinVar(),
                        brca.brcApart1().Pathogenicity_expert(),
                        brca.brcApart1().Allele_frequency_1000_Genomes(),
                        brca.brcApart1().Functional_analysis_result_LOVD(),
                        brca.brcApart1().AMR_Allele_frequency_1000_Genomes(),
                        brca.brcApart1().Variant_in_ESP(),
                        brca.brcApart1().Variant_in_BIC(),
                        brca.brcApart1().Clinical_significance_ENIGMA(),
                        brca.brcApart1().Max_Allele_Frequency(),
                        brca.brcApart1().Allele_count_AMR_ExAC(),
                        brca.brcApart1().Variant_in_ENIGMA(),
                        brca.brcApart1().BX_ID_ESP(),
                        brca.brcApart1().Patient_nationality_BIC(),
                        brca.brcApart1().BX_ID_1000_Genomes(),
                        brca.brcApart1().Genomic_Coordinate_hg37(),
                        brca.brcApart1().Genomic_Coordinate_hg36(),
                        brca.brcApart1().EUR_Allele_frequency_1000_Genomes(),
                        brca.brcApart1().Number_of_family_member_carrying_mutation_BIC(),
                        brca.brcApart1().Segregation_LR_exLOVD(),
                        brca.brcApart1().Allele_Frequency(),
                        brca.brcApart1().Minor_allele_frequency_percent_ESP(),
                        brca.brcApart1().Allele_frequency_ExAC(),
                        brca.brcApart1().Mutation_type_BIC(),
                        brca.brcApart1().Assertion_method_citation_ENIGMA(),
                        brca.brcApart1().Condition_ID_type_ENIGMA(),
                        brca.brcApart1().Allele_count_OTH_ExAC(),
                        brca.brcApart1().HGVS_protein_LOVD(),
                        brca.brcApart1().Variant_in_ClinVar(),
                        brca.brcApart1().Clinical_importance_BIC(),
                        brca.brcApart1().Discordant(),
                        brca.brcApart2().Allele_count_FIN_ExAC(),
                        brca.brcApart2().Condition_category_ENIGMA(),
                        brca.brcApart2().Allele_Frequency_ESP(),
                        brca.brcApart2().Homozygous_count_OTH_ExAC(),
                        brca.brcApart2().Genetic_origin_LOVD(),
                        brca.brcApart2().Homozygous_count_AMR_ExAC(),
                        brca.brcApart2().Clinical_Significance_ClinVar(),
                        brca.brcApart2().AA_Allele_Frequency_ESP(),
                        brca.brcApart2().Protein_Change(),
                        brca.brcApart2().Variant_in_exLOVD(),
                        brca.brcApart2().EA_Allele_Frequency_ESP(),
                        brca.brcApart2().HGVS_RNA(),
                        brca.brcApart2().Clinical_significance_citations_ENIGMA(),
                        brca.brcApart2().Variant_effect_LOVD(),
                        brca.brcApart2().Polyphen_Prediction(),
                        brca.brcApart2().Data_Release_id(),
                        brca.brcApart2().Hg37_Start(),
                        brca.brcApart2().Hg36_Start(),
                        brca.brcApart2().Sift_Score(),
                        brca.brcApart2().Genomic_Coordinate_hg38(),
                        brca.brcApart2().Alt(),
                        brca.brcApart2().Literature_citation_BIC(),
                        brca.brcApart2().Variant_haplotype_LOVD(),
                        brca.brcApart2().Allele_frequency_NFE_ExAC(),
                        brca.brcApart2().Hg38_Start(),
                        brca.brcApart2().Pos(),
                        brca.brcApart2().Date_last_evaluated_ENIGMA(),
                        brca.brcApart2().Allele_number_SAS_ExAC(),
                        brca.brcApart2().Allele_number_AMR_ExAC(),
                        viccEntryId)
                .execute();
    }

    private void importCGIinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Cgi cgi = (Cgi) object;
        int id = context.insertInto(CGI,
                CGI.TARGETING,
                CGI.SOURCE,
                CGI.PRIMARYTUMORTYPE,
                CGI.DRUGSFULLNAME,
                CGI.CURATOR,
                CGI.DRUGFAMILY,
                CGI.ALTERATION,
                CGI.DRUG,
                CGI.BIOMARKER,
                CGI.DRUGSTATUS,
                CGI.GENE,
                CGI.ASSAYTYPE,
                CGI.ALTERATIONTYPE,
                CGI.EVIDENCELEVEL,
                CGI.ASSOCIATION,
                CGI.METASTATICTUMORTYPE,
                CGI.VICCENTRYID)
                .values(cgi.targeting(),
                        cgi.source(),
                        cgi.primary_tumor_type(),
                        cgi.drugsFullName(),
                        cgi.curator(),
                        cgi.drug_family(),
                        cgi.alteration(),
                        cgi.drug(),
                        cgi.biomarker(),
                        cgi.drug_status(),
                        cgi.gene(),
                        cgi.assay_type(),
                        cgi.alteration_type(),
                        cgi.evidence_level(),
                        cgi.association(),
                        cgi.metastatic_Tumor_Type(),
                        viccEntryId)
                .returning(CGI.ID)
                .fetchOne()
                .getValue(CGI.ID);

        for (String cDNA : cgi.cDNA()) {
            context.insertInto(CGICDNA, CGICDNA.CDNA, CGICDNA.CGIID).values(cDNA, id).execute();
        }

        for (String individualMutation : cgi.individual_mutation()) {
            context.insertInto(CGIINDIVIDUALMUTATION, CGIINDIVIDUALMUTATION.INDIVIDUALMUTATION, CGIINDIVIDUALMUTATION.CGIID)
                    .values(individualMutation, id)
                    .execute();
        }

        for (String gDNA : cgi.gDNA()) {
            context.insertInto(CGIGDNA, CGIGDNA.GDNA, CGIGDNA.CGIID).values(gDNA, id).execute();
        }

        for (String transcript : cgi.transcript()) {
            context.insertInto(CGITRANSCRIPT, CGITRANSCRIPT.TRANSCRIPT, CGITRANSCRIPT.CGIID).values(transcript, id).execute();
        }

        for (String strand : cgi.strand()) {
            context.insertInto(CGISTRAND, CGISTRAND.STRAND, CGISTRAND.CGIID).values(strand, id).execute();
        }

        for (String info : cgi.info()) {
            context.insertInto(CGIINFO, CGIINFO.INFO, CGIINFO.CGIID).values(info, id).execute();
        }

        for (String region : cgi.region()) {
            context.insertInto(CGIREGION, CGIREGION.REGION, CGIREGION.CGIID).values(region, id).execute();
        }
    }

    private void importJAXinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Jax jax = (Jax) object;

        int id = context.insertInto(JAX,
                JAX.RESPONSETYPE,
                JAX.APPROVALSTATUS,
                JAX.EVIDENCETYPE,
                JAX.EFFICACYEVIDENCE,
                JAX.IDJAXSOURCE,
                JAX.VICCENTRYID)
                .values(jax.responseType(), jax.approvalStatus(), jax.evidenceType(), jax.efficacyEvidence(), jax.id(), viccEntryId)
                .returning(JAX.ID)
                .fetchOne()
                .getValue(JAX.ID);

        context.insertInto(JAXMOLECULARPROFILE,
                JAXMOLECULARPROFILE.PROFILENAME,
                JAXMOLECULARPROFILE.IDMOLECULARPROFILE,
                JAXMOLECULARPROFILE.JAXID).values(jax.molecularProfile().profileName(), jax.molecularProfile().id(), id).execute();

        context.insertInto(JAXTHERAPY, JAXTHERAPY.THERAPYNAME, JAXTHERAPY.IDTHERAPY, JAXTHERAPY.JAXID)
                .values(jax.therapy().therapyName(), jax.therapy().id(), id)
                .execute();

        context.insertInto(JAXINDICATIONS, JAXINDICATIONS.SOURCE, JAXINDICATIONS.IDINDICATIONS, JAXINDICATIONS.NAME, JAXINDICATIONS.JAXID)
                .values(jax.indications().source(), jax.indications().id(), jax.indications().name(), id)
                .execute();

        for (JaxReferences references : jax.references()) {
            context.insertInto(JAXREFERENCES,
                    JAXREFERENCES.URL,
                    JAXREFERENCES.IDREFERENCES,
                    JAXREFERENCES.PUBMEDID,
                    JAXREFERENCES.TITLE,
                    JAXREFERENCES.JAXID).values(references.url(), references.id(), references.pubMedId(), references.title(), id).execute();
        }
    }

    private void importJaxTrialsinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        JaxTrials jaxTrials = (JaxTrials) object;

        int id = context.insertInto(JAXTRIALS,
                JAXTRIALS.TITLE,
                JAXTRIALS.GENDER,
                JAXTRIALS.NCTID,
                JAXTRIALS.SPONSORS,
                JAXTRIALS.RECRUITMENT,
                JAXTRIALS.VARIANTREQUIREMENTS,
                JAXTRIALS.UPDATEDATE,
                JAXTRIALS.PHASE,
                JAXTRIALS.VICCENTRYID)
                .values(jaxTrials.title(),
                        jaxTrials.gender(),
                        jaxTrials.nctId(),
                        jaxTrials.sponsors(),
                        jaxTrials.recruitment(),
                        jaxTrials.variantRequirements(),
                        jaxTrials.updateDate(),
                        jaxTrials.updateDate(),
                        viccEntryId)
                .returning(JAXTRIALS.ID)
                .fetchOne()
                .getValue(JAXTRIALS.ID);

        for (JaxTrialsIndications indications : jaxTrials.indications()) {
            context.insertInto(JAXTRIALSINDICATIONS,
                    JAXTRIALSINDICATIONS.SOURCE,
                    JAXTRIALSINDICATIONS.IDINDICATIONS,
                    JAXTRIALSINDICATIONS.NAME,
                    JAXTRIALSINDICATIONS.JAXTRIALSID).values(indications.source(), indications.id(), indications.name(), id).execute();
        }

        for (JaxTrialsVariantRequirementDetails variantRequirementDetails : jaxTrials.variantRequirementDetails()) {
            int id1 = context.insertInto(JAXTRIALSVARIANTREQUIREMENTDETAILS,
                    JAXTRIALSVARIANTREQUIREMENTDETAILS.REQUIREMENTTYPE,
                    JAXTRIALSVARIANTREQUIREMENTDETAILS.JAXTRIALSID)
                    .values(variantRequirementDetails.requirementType(), id)
                    .returning(JAXTRIALSVARIANTREQUIREMENTDETAILS.ID)
                    .fetchOne()
                    .getValue(JAXTRIALSVARIANTREQUIREMENTDETAILS.ID);

            for (JaxTrialsMolecularProfile molecularProfile : variantRequirementDetails.molecularProfiles()) {
                context.insertInto(JAXTRIALSMOLECULARPROFILE,
                        JAXTRIALSMOLECULARPROFILE.PROFILENAME,
                        JAXTRIALSMOLECULARPROFILE.IDMOLECULARPROFILE,
                        JAXTRIALSMOLECULARPROFILE.JAXTRIALSVARIANTREQUIREMENTDETAILSID)
                        .values(molecularProfile.profileName(), molecularProfile.id(), id1)
                        .execute();
            }
        }

        for (JaxTrialsTherapies therapies : jaxTrials.therapies()) {
            context.insertInto(JAXTRIALSTHERAPIES,
                    JAXTRIALSTHERAPIES.IDTHERAPIES,
                    JAXTRIALSTHERAPIES.THERAPYNAME,
                    JAXTRIALSTHERAPIES.JAXTRIALSID).values(therapies.id(), therapies.therapyName(), id).execute();
        }

    }

    private void importPmkbinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Pmkb pmkb = (Pmkb) object;

        int id = context.insertInto(PMKB, PMKB.VICCENTRYID).values(viccEntryId).returning(PMKB.ID).fetchOne().getValue(PMKB.ID);

        for (PmkbTumor tumor : pmkb.tumor()) {
            context.insertInto(PMKBTUMOR, PMKBTUMOR.IDTUMOR, PMKBTUMOR.NAME, PMKBTUMOR.VICCENTRYID)
                    .values(tumor.id(), tumor.name(), id)
                    .execute();
        }

        for (PmkbTissue tissue : pmkb.tissue()) {
            context.insertInto(PMKBTISSUE, PMKBTISSUE.IDTISSUE, PMKBTISSUE.NAME, PMKBTISSUE.VICCENTRYID)
                    .values(tissue.id(), tissue.name(), id)
                    .execute();
        }

        for (PmkbVariant variant : pmkb.variant()) {
            int idVariant = context.insertInto(PMKBVARIANT,
                    PMKBVARIANT.AMINOACIDCHANGE,
                    PMKBVARIANT.GERMLINE,
                    PMKBVARIANT.PARTNERGENE,
                    PMKBVARIANT.CODONS,
                    PMKBVARIANT.DESCRIPTION,
                    PMKBVARIANT.EXONS,
                    PMKBVARIANT.NOTES,
                    PMKBVARIANT.COSMIC,
                    PMKBVARIANT.EFFECT,
                    PMKBVARIANT.CNVTYPE,
                    PMKBVARIANT.IDVARIANT,
                    PMKBVARIANT.CYTOBAND,
                    PMKBVARIANT.VARIANTTYPE,
                    PMKBVARIANT.DNACHANGE,
                    PMKBVARIANT.COORDINATES,
                    PMKBVARIANT.CHROMOSOMEBASEDCNV,
                    PMKBVARIANT.TRANSCRIPT,
                    PMKBVARIANT.DESCRIPTIONTYPE,
                    PMKBVARIANT.CHROMOSOME,
                    PMKBVARIANT.NAME,
                    PMKBVARIANT.VICCENTRYID)
                    .values(variant.aminoAcidChange(),
                            variant.germline(),
                            variant.partnerGene(),
                            variant.codons(),
                            variant.description(),
                            variant.exons(),
                            variant.notes(),
                            variant.cosmic(),
                            variant.effect(),
                            variant.cnvType(),
                            variant.id(),
                            variant.cytoband(),
                            variant.variantType(),
                            variant.dnaChange(),
                            variant.coordinates(),
                            variant.chromosomeBasedCnv(),
                            variant.transcript(),
                            variant.descriptionType(),
                            variant.chromosome(),
                            variant.name(),
                            viccEntryId)
                    .returning(PMKBVARIANT.ID)
                    .fetchOne()
                    .getValue(PMKBVARIANT.ID);

            for (PmkbGene gene : variant.gene()) {
                context.insertInto(PMKBGENE,
                        PMKBGENE.DESCRIPTION,
                        PMKBGENE.CREATEDAT,
                        PMKBGENE.UPDATEDAT,
                        PMKBGENE.ACTIVEIND,
                        PMKBGENE.EXTERNALID,
                        PMKBGENE.IDGENE,
                        PMKBGENE.NAME,
                        PMKBGENE.PMKBVARIANTID)
                        .values(gene.description(),
                                gene.createdAt(),
                                gene.updatedAt(),
                                gene.activeInd(),
                                gene.externalId(),
                                gene.id(),
                                gene.name(),
                                idVariant)
                        .execute();
            }
        }
    }

    private void importOncokbinSQL(int viccEntryId, @NotNull KbSpecificObject object) {
        Oncokb oncokb = (Oncokb) object;
        int id = context.insertInto(ONCOKB, ONCOKB.VICCENTRYID).values(viccEntryId).returning(ONCOKB.ID).fetchOne().getValue(ONCOKB.ID);

        OncoKbBiological oncokbBiological = oncokb.oncoKbBiological();
        if (oncokbBiological != null) {
            int idBiological = context.insertInto(ONCOKBBIOLOGICAL,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECTPMIDS,
                    ONCOKBBIOLOGICAL.ISOFORM,
                    ONCOKBBIOLOGICAL.ENTREZGENEID,
                    ONCOKBBIOLOGICAL.ONCOGENIC,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECT,
                    ONCOKBBIOLOGICAL.REFSEQ,
                    ONCOKBBIOLOGICAL.GENE,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECTABSTRACTS,
                    ONCOKBBIOLOGICAL.VICCENTRYID)
                    .values(oncokbBiological.mutationEffectPmids(),
                            oncokbBiological.Isoform(),
                            oncokbBiological.entrezGeneID(),
                            oncokbBiological.oncogenic(),
                            oncokbBiological.mutationEffect(),
                            oncokbBiological.RefSeq(),
                            oncokbBiological.gene(),
                            oncokbBiological.mutationEffectAbstracts(),
                            id)
                    .returning(ONCOKBBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBBIOLOGICAL.ID);

            OncokbVariant oncokbVariant = oncokb.oncoKbBiological().oncokbVariant();

            int idVariant = context.insertInto(ONCOKBVARIANTBIOLOGICAL,
                    ONCOKBVARIANTBIOLOGICAL.VARIANTRESIDUES,
                    ONCOKBVARIANTBIOLOGICAL.PROTEINSTART,
                    ONCOKBVARIANTBIOLOGICAL.NAME,
                    ONCOKBVARIANTBIOLOGICAL.PROTEINEND,
                    ONCOKBVARIANTBIOLOGICAL.REFRESIDUES,
                    ONCOKBVARIANTBIOLOGICAL.ALTERATION,
                    ONCOKBVARIANTBIOLOGICAL.ONCOKBBIOLOGICALID)
                    .values(oncokbVariant.variantResidues(),
                            oncokbVariant.proteinStart(),
                            oncokbVariant.name(),
                            oncokbVariant.proteinEnd(),
                            oncokbVariant.refResidues(),
                            oncokbVariant.alteration(),
                            idBiological)
                    .returning(ONCOKBVARIANTBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBVARIANTBIOLOGICAL.ID);

            OncoKbConsequence oncoKbConsequence = oncokb.oncoKbBiological().oncokbVariant().oncoKbConsequence();

            context.insertInto(ONCOKBCONSEQUENCESBIOLOGICAL,
                    ONCOKBCONSEQUENCESBIOLOGICAL.TERM,
                    ONCOKBCONSEQUENCESBIOLOGICAL.DESCRIPTION,
                    ONCOKBCONSEQUENCESBIOLOGICAL.ISGENERALLYTRUNCATING,
                    ONCOKBCONSEQUENCESBIOLOGICAL.ONCOKBVARIANTBIOLOGICALID)
                    .values(oncoKbConsequence.term(), oncoKbConsequence.description(), oncoKbConsequence.isGenerallyTruncating(), idVariant)
                    .execute();

            OncokbGene oncoKbGene = oncokb.oncoKbBiological().oncokbVariant().oncokbGene();

            int idGene = context.insertInto(ONCOKBGENEBIOLOGICAL,
                    ONCOKBGENEBIOLOGICAL.ONCOGENE,
                    ONCOKBGENEBIOLOGICAL.NAME,
                    ONCOKBGENEBIOLOGICAL.HUGOSYMBOL,
                    ONCOKBGENEBIOLOGICAL.CURATEDREFSEQ,
                    ONCOKBGENEBIOLOGICAL.ENTREZGENEID,
                    ONCOKBGENEBIOLOGICAL.TSG,
                    ONCOKBGENEBIOLOGICAL.CURATEDISOFORM,
                    ONCOKBGENEBIOLOGICAL.ONCOKBBIOLOGICALID)
                    .values(oncoKbGene.oncogene(),
                            oncoKbGene.name(),
                            oncoKbGene.hugoSymbol(),
                            oncoKbGene.curatedRefSeq(),
                            oncoKbGene.entrezGeneId(),
                            oncoKbGene.tsg(),
                            oncoKbGene.curatedIsoform(),
                            idVariant)
                    .returning(ONCOKBGENEBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBGENEBIOLOGICAL.ID);

            for (String geneAliases : oncoKbGene.geneAliases()) {
                context.insertInto(ONCOKBGENEALIASESBIOLOGICAL,
                        ONCOKBGENEALIASESBIOLOGICAL.GENEALIASES,
                        ONCOKBGENEALIASESBIOLOGICAL.ONCOKBGENEBIOLOGICALID).values(geneAliases, idGene);
            }
        }

        OncoKbClinical oncokbClinical = oncokb.oncoKbClinical();
        if (oncokbClinical != null) {

            int idClinical = context.insertInto(ONCOKBCLINICAL,
                    ONCOKBCLINICAL.REFSEQ,
                    ONCOKBCLINICAL.LEVEL,
                    ONCOKBCLINICAL.ENTREZGENEID,
                    ONCOKBCLINICAL.DRUGPMIDS,
                    ONCOKBCLINICAL.CANCERTYPE,
                    ONCOKBCLINICAL.DRUG,
                    ONCOKBCLINICAL.GENE,
                    ONCOKBCLINICAL.LEVELLABEL,
                    ONCOKBCLINICAL.VICCENTRYID)
                    .values(oncokbClinical.RefSeq(),
                            oncokbClinical.level(),
                            oncokbClinical.entrezGeneID(),
                            oncokbClinical.drugPmids(),
                            oncokbClinical.cancerType(),
                            oncokbClinical.drug(),
                            oncokbClinical.gene(),
                            oncokbClinical.levelLabel(),
                            id)
                    .returning(ONCOKBCLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBCLINICAL.ID);

            for (OncoKbDrugAbstracts drugAbstracts : oncokbClinical.oncoKbDrugAbstracts()) {
                context.insertInto(ONCOKBDRUGABSTRACTSCLINICAL,
                        ONCOKBDRUGABSTRACTSCLINICAL.TEXT,
                        ONCOKBDRUGABSTRACTSCLINICAL.LINK,
                        ONCOKBDRUGABSTRACTSCLINICAL.ONCOKBCLINICALID)
                        .values(drugAbstracts.text(), drugAbstracts.link(), idClinical)
                        .execute();
            }

            OncokbVariant variantClinical = oncokbClinical.oncokbVariant();
            int idClinicalVariant = context.insertInto(ONCOKBVARIANTCLINICAL,
                    ONCOKBVARIANTCLINICAL.VARIANTRESIDUES,
                    ONCOKBVARIANTCLINICAL.PROTEINSTART,
                    ONCOKBVARIANTCLINICAL.NAME,
                    ONCOKBVARIANTCLINICAL.PROTEINEND,
                    ONCOKBVARIANTCLINICAL.REFRESIDUES,
                    ONCOKBVARIANTCLINICAL.ALTERATION,
                    ONCOKBVARIANTCLINICAL.ONCOKBCLINICALID)
                    .values(variantClinical.variantResidues(),
                            variantClinical.proteinStart(),
                            variantClinical.name(),
                            variantClinical.proteinEnd(),
                            variantClinical.refResidues(),
                            variantClinical.alteration(),
                            idClinical)
                    .returning(ONCOKBVARIANTCLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBVARIANTCLINICAL.ID);

            OncoKbConsequence consequenceaClinical = oncokbClinical.oncokbVariant().oncoKbConsequence();

            context.insertInto(ONCOKBCONSEQUENCESCLINICAL,
                    ONCOKBCONSEQUENCESCLINICAL.TERM,
                    ONCOKBCONSEQUENCESCLINICAL.DESCRIPTION,
                    ONCOKBCONSEQUENCESCLINICAL.ISGENERALLYTRUNCATING,
                    ONCOKBCONSEQUENCESCLINICAL.ONCOKBVARIANTCLINICALID)
                    .values(consequenceaClinical.term(),
                            consequenceaClinical.description(),
                            consequenceaClinical.isGenerallyTruncating(),
                            idClinicalVariant)
                    .execute();

            OncokbGene geneClinical = oncokbClinical.oncokbVariant().oncokbGene();

            int idGeneClinical = context.insertInto(ONCOKBGENECLINICAL,
                    ONCOKBGENECLINICAL.ONCOGENE,
                    ONCOKBGENECLINICAL.NAME,
                    ONCOKBGENECLINICAL.HUGOSYMBOL,
                    ONCOKBGENECLINICAL.CURATEDREFSEQ,
                    ONCOKBGENECLINICAL.ENTREZGENEID,
                    ONCOKBGENECLINICAL.TSG,
                    ONCOKBGENECLINICAL.CURATEDISOFORM,
                    ONCOKBGENECLINICAL.ONCOKBCLINICALID)
                    .values(geneClinical.oncogene(),
                            geneClinical.name(),
                            geneClinical.hugoSymbol(),
                            geneClinical.curatedRefSeq(),
                            geneClinical.entrezGeneId(),
                            geneClinical.tsg(),
                            geneClinical.curatedIsoform(),
                            idClinicalVariant)
                    .returning(ONCOKBGENECLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBGENECLINICAL.ID);

            for (String geneAliasesClinical : geneClinical.geneAliases()) {
                context.insertInto(ONCOKBGENEALIASESCLINICAL,
                        ONCOKBGENEALIASESCLINICAL.GENEALIASES,
                        ONCOKBGENEALIASESCLINICAL.ONCOKBGENECLINICALID).values(geneAliasesClinical, idGeneClinical).execute();
            }

        }

    }

    private void importMolecularMatchTrials(int viccEntryId, @NotNull KbSpecificObject object) {
        MolecularMatchTrials molecularMatchTrials = (MolecularMatchTrials) object;
        int id = context.insertInto(MOLECULARMATCHTRIALS,
                MOLECULARMATCHTRIALS.STATUS,
                MOLECULARMATCHTRIALS.STARTDATE,
                MOLECULARMATCHTRIALS.TITLE,
                MOLECULARMATCHTRIALS.SCORE,
                MOLECULARMATCHTRIALS.BRIEFTITLE,
                MOLECULARMATCHTRIALS.LINK,
                MOLECULARMATCHTRIALS.PHASE,
                MOLECULARMATCHTRIALS.IDMOLECULARMATCHTRIALS,
                MOLECULARMATCHTRIALS.STUDYTYPE,
                MOLECULARMATCHTRIALS.VICCENTRYID)
                .values(molecularMatchTrials.status(),
                        molecularMatchTrials.startDate(),
                        molecularMatchTrials.title(),
                        molecularMatchTrials.score(),
                        molecularMatchTrials.briefTitle(),
                        molecularMatchTrials.link(),
                        molecularMatchTrials.phase(),
                        molecularMatchTrials.id(),
                        molecularMatchTrials.studyType(),
                        viccEntryId)
                .returning(MOLECULARMATCHTRIALS.ID)
                .fetchOne()
                .getValue(MOLECULARMATCHTRIALS.ID);

        for (String molecularAlterations : molecularMatchTrials.molecularAlterations()) {
            context.insertInto(MOLECULARMATCHTRIALSALTERATIONS,
                    MOLECULARMATCHTRIALSALTERATIONS.MOLECULARALTERATIONS,
                    MOLECULARMATCHTRIALSALTERATIONS.MOLECULARMATCHTRIALSID).values(molecularAlterations, id).execute();
        }

        for (MolecularMatchTrialsIntervation intervation : molecularMatchTrials.intervation()) {
            int idIntervation = context.insertInto(MOLECULARMATCHTRIALSINTERVATIONS,
                    MOLECULARMATCHTRIALSINTERVATIONS.INTERVENTION_NAME,
                    MOLECULARMATCHTRIALSINTERVATIONS.DESCRIPTION,
                    MOLECULARMATCHTRIALSINTERVATIONS.INTERVENTION_TYPE,
                    MOLECULARMATCHTRIALSINTERVATIONS.MOLECULARMATCHTRIALSID)
                    .values(intervation.intervention_name(), intervation.description(), intervation.intervention_type(), viccEntryId)
                    .returning(MOLECULARMATCHTRIALSINTERVATIONS.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHTRIALSINTERVATIONS.ID);

            if (intervation.other_name() != null) {
                for (String otherName : intervation.other_name()) {
                    context.insertInto(MOLECULARMATCHTRIALSOTHERNAME,
                            MOLECULARMATCHTRIALSOTHERNAME.OTHER_NAME,
                            MOLECULARMATCHTRIALSOTHERNAME.MOLECULARMATCHTRIALSINTERVATIONSID).values(otherName, idIntervation).execute();
                }
            }

            if (intervation.arm_group_label() != null) {
                for (String armGroupLabel : intervation.arm_group_label()) {
                    context.insertInto(MOLECULARMATCHTRIALSOTHERGROUPLABEL,
                            MOLECULARMATCHTRIALSOTHERGROUPLABEL.ARM_GROUP_LABEL,
                            MOLECULARMATCHTRIALSOTHERGROUPLABEL.MOLECULARMATCHTRIALSINTERVATIONSID)
                            .values(armGroupLabel, idIntervation)
                            .execute();
                }
            }
        }

        for (MolecularMatchTrialsLocations locations : molecularMatchTrials.locations()) {
            int idLocations = context.insertInto(MOLECULARMATCHTRIALSLOCATIONS,
                    MOLECULARMATCHTRIALSLOCATIONS.STATUS,
                    MOLECULARMATCHTRIALSLOCATIONS.LAST_NAME,
                    MOLECULARMATCHTRIALSLOCATIONS.EMAIL,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.EMAIL_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.LAST_NAME_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_EXT_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_EXT,
                    MOLECULARMATCHTRIALSLOCATIONS.CITY,
                    MOLECULARMATCHTRIALSLOCATIONS.VALID,
                    MOLECULARMATCHTRIALSLOCATIONS.ZIP,
                    MOLECULARMATCHTRIALSLOCATIONS.CREATED,
                    MOLECULARMATCHTRIALSLOCATIONS.COUNTRY,
                    MOLECULARMATCHTRIALSLOCATIONS.NUMBER,
                    MOLECULARMATCHTRIALSLOCATIONS.IDLOCATIONS,
                    MOLECULARMATCHTRIALSLOCATIONS.LASTUPDATED,
                    MOLECULARMATCHTRIALSLOCATIONS.STATE,
                    MOLECULARMATCHTRIALSLOCATIONS.STREET,
                    MOLECULARMATCHTRIALSLOCATIONS.PO_BOX,
                    MOLECULARMATCHTRIALSLOCATIONS.FAILEDGEOCODE,
                    MOLECULARMATCHTRIALSLOCATIONS.VALIDMESSAGE,
                    MOLECULARMATCHTRIALSLOCATIONS.NAME,
                    MOLECULARMATCHTRIALSLOCATIONS.MOLECULARMATCHTRIALSID)
                    .values(locations.status(),
                            locations.last_name(),
                            locations.email(),
                            locations.phone(),
                            locations.phone_backup(),
                            locations.email_backup(),
                            locations.last_name_backup(),
                            locations.phone_ext_backup(),
                            locations.phone_ext(),
                            locations.city(),
                            locations.valid(),
                            locations.zip(),
                            locations.created(),
                            locations.country(),
                            locations.number(),
                            locations.id(),
                            locations.lastUpdated(),
                            locations.state(),
                            locations.street(),
                            locations.po_box(),
                            locations.failedGeocode(),
                            locations.validMessage(),
                            locations.name(),
                            id)
                    .returning(MOLECULARMATCHTRIALSLOCATIONS.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHTRIALSLOCATIONS.ID);

            if (locations.contact() != null) {
                context.insertInto(MOLECULARMATCHTRIALSCONTACT,
                        MOLECULARMATCHTRIALSCONTACT.PHONE,
                        MOLECULARMATCHTRIALSCONTACT.NAME,
                        MOLECULARMATCHTRIALSCONTACT.EMAIL,
                        MOLECULARMATCHTRIALSCONTACT.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.contact().phone(), locations.contact().name(), locations.contact().email(), idLocations)
                        .execute();
            }

            if (locations.location() != null) {
                int idLocation = context.insertInto(MOLECULARMATCHTRIALSLOCATION,
                        MOLECULARMATCHTRIALSLOCATION.TYPE,
                        MOLECULARMATCHTRIALSLOCATION.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.location().type(), idLocations)
                        .returning(MOLECULARMATCHTRIALSLOCATION.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHTRIALSLOCATION.ID);

                for (String coordinates : locations.location().coordinates()) {
                    context.insertInto(MOLECULARMATCHTRIALSCOORDINATES,
                            MOLECULARMATCHTRIALSCOORDINATES.COORDINATES,
                            MOLECULARMATCHTRIALSCOORDINATES.MOLECULARMATCHTRIALSLOCATIONID).values(coordinates, idLocation).execute();
                }
            }

            if (locations.geo() != null) {
                context.insertInto(MOLECULARMATCHTRIALSGEO,
                        MOLECULARMATCHTRIALSGEO.LAT,
                        MOLECULARMATCHTRIALSGEO.LON,
                        MOLECULARMATCHTRIALSGEO.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.geo().lat(), locations.geo().lon(), idLocations)
                        .execute();
            }
        }

        if (molecularMatchTrials.overallContact() != null) {
            context.insertInto(MOLECULARMATCHTRIALSOVERALLCONTACT,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.PHONE,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.LAST_NAME,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.PHONE_EXT,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.COUNTRY,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.EMAIL,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.AFFILIATION,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.CITY,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.NAME,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.ZIP,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.URL,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.STREET,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.TYPE,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.MOLECULARMATCHTRIALSID)
                    .values(molecularMatchTrials.overallContact().phone(),
                            molecularMatchTrials.overallContact().last_name(),
                            molecularMatchTrials.overallContact().phone_ext(),
                            molecularMatchTrials.overallContact().country(),
                            molecularMatchTrials.overallContact().email(),
                            molecularMatchTrials.overallContact().affiliation(),
                            molecularMatchTrials.overallContact().city(),
                            molecularMatchTrials.overallContact().name(),
                            molecularMatchTrials.overallContact().zip(),
                            molecularMatchTrials.overallContact().url(),
                            molecularMatchTrials.overallContact().street(),
                            molecularMatchTrials.overallContact().type(),
                            id)
                    .execute();
        }

        for (MolecularMatchTrialsTags tags : molecularMatchTrials.tags()) {
            context.insertInto(MOLECULARMATCHTRIALSTAGS,
                    MOLECULARMATCHTRIALSTAGS.FACET,
                    MOLECULARMATCHTRIALSTAGS.COMPOSITEKEY,
                    MOLECULARMATCHTRIALSTAGS.SUPPRESS,
                    MOLECULARMATCHTRIALSTAGS.FILTERTYPE,
                    MOLECULARMATCHTRIALSTAGS.TERM,
                    MOLECULARMATCHTRIALSTAGS.CUSTOM,
                    MOLECULARMATCHTRIALSTAGS.PRIORITY,
                    MOLECULARMATCHTRIALSTAGS.ALIAS,
                    MOLECULARMATCHTRIALSTAGS.MANUALSUPPRESS,
                    MOLECULARMATCHTRIALSTAGS.GENERATEDBY,
                    MOLECULARMATCHTRIALSTAGS.GENERATEDBYTERM,
                    MOLECULARMATCHTRIALSTAGS.IDTAGS,
                    MOLECULARMATCHTRIALSTAGS.MANUALPRIORITY,
                    MOLECULARMATCHTRIALSTAGS.MOLECULARMATCHTRIALSID)
                    .values(tags.facet(),
                            tags.compositeKey(),
                            tags.suppress(),
                            tags.filterType(),
                            tags.term(),
                            tags.custom(),
                            tags.priority(),
                            tags.alias(),
                            tags.manualSuppress(),
                            tags.generatedBy(),
                            tags.generatedByTerm(),
                            tags.id(),
                            tags.manualPriority(),
                            id)
                    .execute();
        }

    }

    private void importCivic(int viccEntryId, @NotNull KbSpecificObject object) {
        Civic civic = (Civic) object;

        int id = context.insertInto(CIVIC,
                CIVIC.ENTREZNAME,
                CIVIC.CIVICACTIONABILITYSCORE,
                CIVIC.ALLELEREGISTRYID,
                CIVIC.GENEID,
                CIVIC.NAME,
                CIVIC.ENTREZID,
                CIVIC.TYPE,
                CIVIC.IDCIVIC,
                CIVIC.DESCRIPTION,
                CIVIC.VICCENTRYID)
                .values(civic.entrezName(),
                        civic.civicActionabilityScore(),
                        civic.alleleRegistryId(),
                        civic.geneId(),
                        civic.name(),
                        civic.entrezId(),
                        civic.type(),
                        civic.id(),
                        civic.description(),
                        viccEntryId)
                .returning(CIVIC.ID)
                .fetchOne()
                .getValue(CIVIC.ID);

        for (String assertions : civic.assertions()) {
            context.insertInto(CIVICASSERTIONS, CIVICASSERTIONS.ASSERTIONS, CIVICASSERTIONS.CIVICID).values(assertions, id).execute();
        }

        for (String hgvsExpression : civic.hgvs_expressions()) {
            context.insertInto(CIVICHGVSEXPRESSIONS, CIVICHGVSEXPRESSIONS.HGVS_EXPRESSIONS, CIVICHGVSEXPRESSIONS.CIVICID)
                    .values(hgvsExpression, id)
                    .execute();
        }

        for (String clinvarEntries : civic.clinvarEntries()) {
            context.insertInto(CIVICCLINVARENTRIES, CIVICCLINVARENTRIES.CLINVARENTRIES, CIVICCLINVARENTRIES.CIVICID)
                    .values(clinvarEntries, id)
                    .execute();
        }

        for (String variantAliases : civic.variantAliases()) {
            context.insertInto(CIVICVARIANTALIASES, CIVICVARIANTALIASES.VARIANTALIASES, CIVICVARIANTALIASES.CIVICID)
                    .values(variantAliases, id)
                    .execute();
        }

        for (CivicVariantTypes variantTypes : civic.variantTypes()) {
            context.insertInto(CIVICVARIANTTYPES,
                    CIVICVARIANTTYPES.DISPLAYNAME,
                    CIVICVARIANTTYPES.DESCRIPTION,
                    CIVICVARIANTTYPES.URL,
                    CIVICVARIANTTYPES.SOID,
                    CIVICVARIANTTYPES.IDVARIANTTYPES,
                    CIVICVARIANTTYPES.NAME,
                    CIVICVARIANTTYPES.CIVICID)
                    .values(variantTypes.displayName(),
                            variantTypes.description(),
                            variantTypes.url(),
                            variantTypes.soId(),
                            variantTypes.id(),
                            variantTypes.name(),
                            id)
                    .execute();
        }

        if (civic.provisional_values() != null) {
            context.insertInto(CIVICDESCRIPTION, CIVICDESCRIPTION.REVISIONID, CIVICDESCRIPTION.VALUE, CIVICDESCRIPTION.CIVICID)
                    .values(civic.provisional_values().revision_id(), civic.provisional_values().value(), id)
                    .execute();
        }

        context.insertInto(CIVICCOORDINATES,
                CIVICCOORDINATES.CHROMOSOME2,
                CIVICCOORDINATES.REFERENCEBASES,
                CIVICCOORDINATES.START2,
                CIVICCOORDINATES.VARIANTBASES,
                CIVICCOORDINATES.STOP,
                CIVICCOORDINATES.STOP2,
                CIVICCOORDINATES.REPRESENTATIVETRANSCRIPT2,
                CIVICCOORDINATES.START,
                CIVICCOORDINATES.REPRESENTATIVETRANSCRIPT,
                CIVICCOORDINATES.ENSEMBLVERSION,
                CIVICCOORDINATES.CHROMOSOME,
                CIVICCOORDINATES.REFERENCEBUILD,
                CIVICCOORDINATES.CIVICID)
                .values(civic.coordinates().chromosome2(),
                        civic.coordinates().referenceBases(),
                        civic.coordinates().start2(),
                        civic.coordinates().variantBases(),
                        civic.coordinates().stop(),
                        civic.coordinates().stop2(),
                        civic.coordinates().representativeTranscript2(),
                        civic.coordinates().start(),
                        civic.coordinates().representativeTranscript(),
                        civic.coordinates().ensemblVersion(),
                        civic.coordinates().chromosome(),
                        civic.coordinates().referenceBuild(),
                        id)
                .execute();

        for (CivicVariantGroup variantGroup : civic.variantGroups()) {
            int idVariantGroup = context.insertInto(CIVICVARIANTSGROUPS,
                    CIVICVARIANTSGROUPS.IDVARIANTSGROUPS,
                    CIVICVARIANTSGROUPS.TYPE,
                    CIVICVARIANTSGROUPS.DESCRIPTION,
                    CIVICVARIANTSGROUPS.NAME,
                    CIVICVARIANTSGROUPS.CIVICID)
                    .values(variantGroup.id(), variantGroup.type(), variantGroup.description(), variantGroup.name(), id)
                    .returning(CIVICVARIANTSGROUPS.ID)
                    .fetchOne()
                    .getValue(CIVICVARIANTSGROUPS.ID);

            for (CivicVariants variants : variantGroup.variants()) {
                int idVariantGroupVariants = context.insertInto(CIVICVARIANTSGROUPSVARIANTS,
                        CIVICVARIANTSGROUPSVARIANTS.ENTREZ_NAME,
                        CIVICVARIANTSGROUPSVARIANTS.DESCRIPTION,
                        CIVICVARIANTSGROUPSVARIANTS.CIVIC_ACTIONABILITY_SCORE,
                        CIVICVARIANTSGROUPSVARIANTS.GENE_ID,
                        CIVICVARIANTSGROUPSVARIANTS.ENTREZ_ID,
                        CIVICVARIANTSGROUPSVARIANTS.TYPE,
                        CIVICVARIANTSGROUPSVARIANTS.IDVARIANTS,
                        CIVICVARIANTSGROUPSVARIANTS.NAME,
                        CIVICVARIANTSGROUPSVARIANTS.CIVICVARIANTSGROUPSID)
                        .values(variants.entrez_name(),
                                variants.description(),
                                variants.civic_actionability_score(),
                                variants.gene_id(),
                                variants.entrez_id(),
                                variants.type(),
                                variants.id(),
                                variants.name(),
                                idVariantGroup)
                        .returning(CIVICVARIANTSGROUPSVARIANTS.ID)
                        .fetchOne()
                        .getValue(CIVICVARIANTSGROUPSVARIANTS.ID);

                if (variants.coordinates() != null) {
                    context.insertInto(CIVICVARIANTSGROUPSCOORDINATES,
                            CIVICVARIANTSGROUPSCOORDINATES.CHROMOSOME2,
                            CIVICVARIANTSGROUPSCOORDINATES.REFERENCEBASES,
                            CIVICVARIANTSGROUPSCOORDINATES.START2,
                            CIVICVARIANTSGROUPSCOORDINATES.VARIANTBASES,
                            CIVICVARIANTSGROUPSCOORDINATES.STOP,
                            CIVICVARIANTSGROUPSCOORDINATES.STOP2,
                            CIVICVARIANTSGROUPSCOORDINATES.REPRESENTATIVETRANSCRIPT2,
                            CIVICVARIANTSGROUPSCOORDINATES.START,
                            CIVICVARIANTSGROUPSCOORDINATES.REPRESENTATIVETRANSCRIPT,
                            CIVICVARIANTSGROUPSCOORDINATES.ENSEMBLVERSION,
                            CIVICVARIANTSGROUPSCOORDINATES.CHROMOSOME,
                            CIVICVARIANTSGROUPSCOORDINATES.REFERENCEBUILD,
                            CIVICVARIANTSGROUPSCOORDINATES.CIVICVARIANTSGROUPSVARIANTSID)
                            .values(variants.coordinates().chromosome2(),
                                    variants.coordinates().referenceBases(),
                                    variants.coordinates().start2(),
                                    variants.coordinates().variantBases(),
                                    variants.coordinates().stop(),
                                    variants.coordinates().stop2(),
                                    variants.coordinates().representativeTranscript2(),
                                    variants.coordinates().start(),
                                    variants.coordinates().representativeTranscript(),
                                    variants.coordinates().ensemblVersion(),
                                    variants.coordinates().chromosome(),
                                    variants.coordinates().referenceBuild(),
                                    idVariantGroupVariants)
                            .execute();
                }

                for (CivicVariantTypes variantTypesGroup : variants.variant_types()) {
                    context.insertInto(CIVICVARIANTSGROUPSTYPES,
                            CIVICVARIANTSGROUPSTYPES.DISPLAYNAME,
                            CIVICVARIANTSGROUPSTYPES.DESCRIPTION,
                            CIVICVARIANTSGROUPSTYPES.URL,
                            CIVICVARIANTSGROUPSTYPES.SOID,
                            CIVICVARIANTSGROUPSTYPES.IDVARIANTTYPES,
                            CIVICVARIANTSGROUPSTYPES.NAME,
                            CIVICVARIANTSGROUPSTYPES.CIVICVARIANTSGROUPSVARIANTSID)
                            .values(variantTypesGroup.displayName(),
                                    variantTypesGroup.description(),
                                    variantTypesGroup.url(),
                                    variantTypesGroup.soId(),
                                    variantTypesGroup.id(),
                                    variantTypesGroup.name(),
                                    idVariantGroupVariants)
                            .execute();
                }
            }
        }

        for (CivicEvidenceItems evidenceItems : civic.evidenceItem()) {
            int idEvidenceItems = context.insertInto(CIVICEVIDENCEITEMS,
                    CIVICEVIDENCEITEMS.STATUS,
                    CIVICEVIDENCEITEMS.RATING,
                    CIVICEVIDENCEITEMS.DRUGINTERACTIONTYPE,
                    CIVICEVIDENCEITEMS.DESCRIPTION,
                    CIVICEVIDENCEITEMS.OPENCHANGECOUNT,
                    CIVICEVIDENCEITEMS.EVIDENCETYPE,
                    CIVICEVIDENCEITEMS.VARIANTORIGIN,
                    CIVICEVIDENCEITEMS.EVIDENCEDIRECTION,
                    CIVICEVIDENCEITEMS.VARIANTID,
                    CIVICEVIDENCEITEMS.CLINICALSIGNIFICANCE,
                    CIVICEVIDENCEITEMS.EVIDENCELEVEL,
                    CIVICEVIDENCEITEMS.TYPE,
                    CIVICEVIDENCEITEMS.IDEVIDENCEITEMS,
                    CIVICEVIDENCEITEMS.NAME,
                    CIVICEVIDENCEITEMS.CIVICID)
                    .values(evidenceItems.status(),
                            evidenceItems.rating(),
                            evidenceItems.drugInteractionType(),
                            evidenceItems.description(),
                            evidenceItems.openChangeCount(),
                            evidenceItems.evidenceType(),
                            evidenceItems.variantOrigin(),
                            evidenceItems.evidenceDirection(),
                            evidenceItems.variantId(),
                            evidenceItems.clinicalSignificance(),
                            evidenceItems.evidenceLevel(),
                            evidenceItems.type(),
                            evidenceItems.id(),
                            evidenceItems.name(),
                            id)
                    .returning(CIVICEVIDENCEITEMS.ID)
                    .fetchOne()
                    .getValue(CIVICEVIDENCEITEMS.ID);

            for (CivicDrugs drugs : evidenceItems.drugs()) {
                context.insertInto(CIVICDRUGS, CIVICDRUGS.PUBCHEMID, CIVICDRUGS.IDDRUGS, CIVICDRUGS.NAME, CIVICDRUGS.CIVICEVIDENCEITEMSID)
                        .values(drugs.pubchemId(), drugs.id(), drugs.name(), idEvidenceItems)
                        .execute();
            }
            context.insertInto(CIVICDISEASE,
                    CIVICDISEASE.DOID,
                    CIVICDISEASE.URL,
                    CIVICDISEASE.DISPLAYNAME,
                    CIVICDISEASE.IDDISEASE,
                    CIVICDISEASE.NAME,
                    CIVICDISEASE.CIVICEVIDENCEITEMSID)
                    .values(evidenceItems.disease().doid(),
                            evidenceItems.disease().url(),
                            evidenceItems.disease().displayName(),
                            evidenceItems.disease().id(),
                            evidenceItems.disease().name(),
                            idEvidenceItems)
                    .execute();

            int idEvidenceItemsSource = context.insertInto(CIVICEVIDENCEITEMSSOURCE,
                    CIVICEVIDENCEITEMSSOURCE.STATUS,
                    CIVICEVIDENCEITEMSSOURCE.OPENACCESS,
                    CIVICEVIDENCEITEMSSOURCE.NAME,
                    CIVICEVIDENCEITEMSSOURCE.JOURNAL,
                    CIVICEVIDENCEITEMSSOURCE.CITATION,
                    CIVICEVIDENCEITEMSSOURCE.PMC_ID,
                    CIVICEVIDENCEITEMSSOURCE.FULLJOURNALTITLE,
                    CIVICEVIDENCEITEMSSOURCE.SOURCEURL,
                    CIVICEVIDENCEITEMSSOURCE.PUBMEDID,
                    CIVICEVIDENCEITEMSSOURCE.ISREVIEW,
                    CIVICEVIDENCEITEMSSOURCE.IDSOURCE,
                    CIVICEVIDENCEITEMSSOURCE.CIVICEVIDENCEITEMSID)
                    .values(evidenceItems.source().status(),
                            evidenceItems.source().openAccess(),
                            evidenceItems.source().name(),
                            evidenceItems.source().journal(),
                            evidenceItems.source().citation(),
                            evidenceItems.source().pmc_Id(),
                            evidenceItems.source().fullJournalTitle(),
                            evidenceItems.source().sourceUrl(),
                            evidenceItems.source().pubmedId(),
                            evidenceItems.source().isReview(),
                            evidenceItems.source().id(),
                            idEvidenceItems)
                    .returning(CIVICEVIDENCEITEMSSOURCE.ID)
                    .fetchOne()
                    .getValue(CIVICEVIDENCEITEMSSOURCE.ID);

            context.insertInto(CIVICEVIDENCEITEMSPUBLICATION,
                    CIVICEVIDENCEITEMSPUBLICATION.YEAR,
                    CIVICEVIDENCEITEMSPUBLICATION.DAY,
                    CIVICEVIDENCEITEMSPUBLICATION.MONTH,
                    CIVICEVIDENCEITEMSPUBLICATION.CIVICEVIDENCEITEMSSOURCEID)
                    .values(evidenceItems.source().publicationDate().year(),
                            evidenceItems.source().publicationDate().day(),
                            evidenceItems.source().publicationDate().month(),
                            idEvidenceItemsSource)
                    .execute();

            for (CivicClinicalTrial clinicalTrial : evidenceItems.source().clinicalTrials()) {
                context.insertInto(CIVICEVIDENCEITEMSCLINICALTRIAL,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.NCT_ID,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.DESCRIPTION,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.CLINICAL_TRIAL_URL,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.NAME,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.CIVICEVIDENCEITEMSSOURCEID)
                        .values(clinicalTrial.nct_id(),
                                clinicalTrial.description(),
                                clinicalTrial.clinical_trial_url(),
                                clinicalTrial.name(),
                                idEvidenceItemsSource)
                        .execute();
            }

        }

        if (civic.sources() != null) {
            for (CivicSource source : civic.sources()) {
                int idSource = context.insertInto(CIVICSOURCE,
                        CIVICSOURCE.STATUS,
                        CIVICSOURCE.OPENACCESS,
                        CIVICSOURCE.NAME,
                        CIVICSOURCE.JOURNAL,
                        CIVICSOURCE.CITATION,
                        CIVICSOURCE.PMC_ID,
                        CIVICSOURCE.FULLJOURNALTITLE,
                        CIVICSOURCE.SOURCEURL,
                        CIVICSOURCE.PUBMEDID,
                        CIVICSOURCE.ISREVIEW,
                        CIVICSOURCE.IDSOURCE,
                        CIVICSOURCE.CIVICID)
                        .values(source.status(),
                                source.openAccess(),
                                source.name(),
                                source.journal(),
                                source.citation(),
                                source.pmc_Id(),
                                source.fullJournalTitle(),
                                source.sourceUrl(),
                                source.pubmedId(),
                                source.isReview(),
                                source.id(),
                                id)
                        .returning(CIVICSOURCE.ID)
                        .fetchOne()
                        .getValue(CIVICSOURCE.ID);

                context.insertInto(CIVICPUBLICATION,
                        CIVICPUBLICATION.YEAR,
                        CIVICPUBLICATION.DAY,
                        CIVICPUBLICATION.MONTH,
                        CIVICPUBLICATION.CIVICSOURCEID)
                        .values(source.publicationDate().year(), source.publicationDate().day(), source.publicationDate().month(), idSource)
                        .execute();

                for (CivicClinicalTrial clinicalTrial : source.clinicalTrials()) {
                    context.insertInto(CIVICCLINICALTRIAL,
                            CIVICCLINICALTRIAL.NCT_ID,
                            CIVICCLINICALTRIAL.DESCRIPTION,
                            CIVICCLINICALTRIAL.CLINICAL_TRIAL_URL,
                            CIVICCLINICALTRIAL.NAME,
                            CIVICCLINICALTRIAL.CIVICSOURCEID)
                            .values(clinicalTrial.nct_id(),
                                    clinicalTrial.description(),
                                    clinicalTrial.clinical_trial_url(),
                                    clinicalTrial.name(),
                                    idSource)
                            .execute();
                }
            }
        }

        context.insertInto(CIVICERROR, CIVICERROR.CIVICID).values(id).execute();

        int idLifeActions = context.insertInto(CIVICLIFECYCLEACTIONS, CIVICLIFECYCLEACTIONS.CIVICID)
                .values(id)
                .returning(CIVICLIFECYCLEACTIONS.ID)
                .fetchOne()
                .getValue(CIVICLIFECYCLEACTIONS.ID);

        if (civic.lifecycleActions().lastCommentedOn() != null) {
            int idLastCommendOn =
                    context.insertInto(CIVICLASTCOMMENTEDON, CIVICLASTCOMMENTEDON.TIMESTAMP, CIVICLASTCOMMENTEDON.CIVICLIFECYCLEACTIONSID)
                            .values(civic.lifecycleActions().lastCommentedOn().timestamp(), idLifeActions)
                            .returning(CIVICLASTCOMMENTEDON.ID)
                            .fetchOne()
                            .getValue(CIVICLASTCOMMENTEDON.ID);

            CivicUser userLastCommend = civic.lifecycleActions().lastCommentedOn().user();
            int idLastCommentUser = context.insertInto(CIVICLASTCOMMENTEDONUSER,
                    CIVICLASTCOMMENTEDONUSER.USERNAME,
                    CIVICLASTCOMMENTEDONUSER.AREAOFEXPERTISE,
                    CIVICLASTCOMMENTEDONUSER.TWITTERHANDLE,
                    CIVICLASTCOMMENTEDONUSER.NAME,
                    CIVICLASTCOMMENTEDONUSER.BIO,
                    CIVICLASTCOMMENTEDONUSER.URL,
                    CIVICLASTCOMMENTEDONUSER.CREATEDAT,
                    CIVICLASTCOMMENTEDONUSER.ACCEPTEDLICENSE,
                    CIVICLASTCOMMENTEDONUSER.AFFILIATION,
                    CIVICLASTCOMMENTEDONUSER.AVATARURL,
                    CIVICLASTCOMMENTEDONUSER.ROLE,
                    CIVICLASTCOMMENTEDONUSER.FACEBOOKPROFILE,
                    CIVICLASTCOMMENTEDONUSER.LINKEDINPROFILE,
                    CIVICLASTCOMMENTEDONUSER.ORCID,
                    CIVICLASTCOMMENTEDONUSER.DISPLAYNAME,
                    CIVICLASTCOMMENTEDONUSER.LASTSEENAT,
                    CIVICLASTCOMMENTEDONUSER.FEATUREDEXPERT,
                    CIVICLASTCOMMENTEDONUSER.IDUSER,
                    CIVICLASTCOMMENTEDONUSER.SIGNUPCOMPLETE,
                    CIVICLASTCOMMENTEDONUSER.CIVICLASTCOMMENTEDONID)
                    .values(userLastCommend.username(),
                            userLastCommend.areaOfExpertise(),
                            userLastCommend.twitterHandle(),
                            userLastCommend.name(),
                            userLastCommend.bio(),
                            userLastCommend.url(),
                            userLastCommend.createdAt(),
                            userLastCommend.acceptedLicense(),
                            userLastCommend.affiliation(),
                            userLastCommend.avatarUrl(),
                            userLastCommend.role(),
                            userLastCommend.facebookProfile(),
                            userLastCommend.linkedinProfile(),
                            userLastCommend.orcid(),
                            userLastCommend.displayName(),
                            userLastCommend.lastSeenAt(),
                            userLastCommend.featuredExpert(),
                            userLastCommend.id(),
                            userLastCommend.signupComplete(),
                            idLastCommendOn)
                    .returning(CIVICLASTCOMMENTEDONUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTCOMMENTEDONUSER.ID);

            context.insertInto(CIVICLASTCOMMENTEDONAVATARS,
                    CIVICLASTCOMMENTEDONAVATARS.X32,
                    CIVICLASTCOMMENTEDONAVATARS.X14,
                    CIVICLASTCOMMENTEDONAVATARS.X64,
                    CIVICLASTCOMMENTEDONAVATARS.X128,
                    CIVICLASTCOMMENTEDONAVATARS.CIVICLASTCOMMENTEDONUSERID)
                    .values(userLastCommend.avatars().x32(),
                            userLastCommend.avatars().x14(),
                            userLastCommend.avatars().x64(),
                            userLastCommend.avatars().x128(),
                            idLastCommentUser)
                    .execute();

            int idLastCommentOnOrganization = context.insertInto(CIVICLASTCOMMENTEDONORGANIZATION,
                    CIVICLASTCOMMENTEDONORGANIZATION.URL,
                    CIVICLASTCOMMENTEDONORGANIZATION.IDORGANIZATION,
                    CIVICLASTCOMMENTEDONORGANIZATION.DESCRIPTION,
                    CIVICLASTCOMMENTEDONORGANIZATION.NAME,
                    CIVICLASTCOMMENTEDONORGANIZATION.CIVICLASTCOMMENTEDONUSERID)
                    .values(userLastCommend.organization().url(),
                            userLastCommend.organization().id(),
                            userLastCommend.organization().description(),
                            userLastCommend.organization().name(),
                            idLastCommentUser)
                    .returning(CIVICLASTCOMMENTEDONORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTCOMMENTEDONORGANIZATION.ID);

            context.insertInto(CIVICLASTCOMMENTEDONPROFILEIMAGE,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X32,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X256,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X14,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X64,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X128,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.CIVICLASTCOMMENTEDONORGANIZATIONID)
                    .values(userLastCommend.organization().profileImage().x32(),
                            userLastCommend.organization().profileImage().x256(),
                            userLastCommend.organization().profileImage().x14(),
                            userLastCommend.organization().profileImage().x64(),
                            userLastCommend.organization().profileImage().x128(),
                            idLastCommentOnOrganization)
                    .execute();
        }

        if (civic.lifecycleActions().lastModified() != null) {
            int idLastModiefied =
                    context.insertInto(CIVICLASTMODIFIED, CIVICLASTMODIFIED.TIMESTAMP, CIVICLASTMODIFIED.CIVICLIFECYCLEACTIONSID)
                            .values(civic.lifecycleActions().lastModified().timestamp(), idLifeActions)
                            .returning(CIVICLASTMODIFIED.ID)
                            .fetchOne()
                            .getValue(CIVICLASTMODIFIED.ID);

            CivicUser userModified = civic.lifecycleActions().lastModified().user();
            int idLastModiefiedUser = context.insertInto(CIVICLASTMODIFIEDUSER,
                    CIVICLASTMODIFIEDUSER.USERNAME,
                    CIVICLASTMODIFIEDUSER.AREAOFEXPERTISE,
                    CIVICLASTMODIFIEDUSER.TWITTERHANDLE,
                    CIVICLASTMODIFIEDUSER.NAME,
                    CIVICLASTMODIFIEDUSER.BIO,
                    CIVICLASTMODIFIEDUSER.URL,
                    CIVICLASTMODIFIEDUSER.CREATEDAT,
                    CIVICLASTMODIFIEDUSER.ACCEPTEDLICENSE,
                    CIVICLASTMODIFIEDUSER.AFFILIATION,
                    CIVICLASTMODIFIEDUSER.AVATARURL,
                    CIVICLASTMODIFIEDUSER.ROLE,
                    CIVICLASTMODIFIEDUSER.FACEBOOKPROFILE,
                    CIVICLASTMODIFIEDUSER.LINKEDINPROFILE,
                    CIVICLASTMODIFIEDUSER.ORCID,
                    CIVICLASTMODIFIEDUSER.DISPLAYNAME,
                    CIVICLASTMODIFIEDUSER.LASTSEENAT,
                    CIVICLASTMODIFIEDUSER.FEATUREDEXPERT,
                    CIVICLASTMODIFIEDUSER.IDUSER,
                    CIVICLASTMODIFIEDUSER.SIGNUPCOMPLETE,
                    CIVICLASTMODIFIEDUSER.CIVICLASTMODIFIEDID)
                    .values(userModified.username(),
                            userModified.areaOfExpertise(),
                            userModified.twitterHandle(),
                            userModified.name(),
                            userModified.bio(),
                            userModified.url(),
                            userModified.createdAt(),
                            userModified.acceptedLicense(),
                            userModified.affiliation(),
                            userModified.avatarUrl(),
                            userModified.role(),
                            userModified.facebookProfile(),
                            userModified.linkedinProfile(),
                            userModified.orcid(),
                            userModified.displayName(),
                            userModified.lastSeenAt(),
                            userModified.featuredExpert(),
                            userModified.id(),
                            userModified.signupComplete(),
                            idLastModiefied)
                    .returning(CIVICLASTMODIFIEDUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTMODIFIEDUSER.ID);

            context.insertInto(CIVICLASTMODIFIEDAVATARS,
                    CIVICLASTMODIFIEDAVATARS.X32,
                    CIVICLASTMODIFIEDAVATARS.X14,
                    CIVICLASTMODIFIEDAVATARS.X64,
                    CIVICLASTMODIFIEDAVATARS.X128,
                    CIVICLASTMODIFIEDAVATARS.CIVICLASTMODIFIEDID)
                    .values(userModified.avatars().x32(),
                            userModified.avatars().x14(),
                            userModified.avatars().x64(),
                            userModified.avatars().x128(),
                            idLastModiefiedUser)
                    .execute();

            int idLastModiefiedOrganization = context.insertInto(CIVICLASTMODIFIEDORGANIZATION,
                    CIVICLASTMODIFIEDORGANIZATION.URL,
                    CIVICLASTMODIFIEDORGANIZATION.IDORGANIZATION,
                    CIVICLASTMODIFIEDORGANIZATION.DESCRIPTION,
                    CIVICLASTMODIFIEDORGANIZATION.NAME,
                    CIVICLASTMODIFIEDORGANIZATION.CIVICLASTMODIFIEDUSERID)
                    .values(userModified.organization().url(),
                            userModified.organization().id(),
                            userModified.organization().description(),
                            userModified.organization().name(),
                            idLastModiefiedUser)
                    .returning(CIVICLASTMODIFIEDORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTMODIFIEDORGANIZATION.ID);

            if (userModified.organization().profileImage() != null) {
                context.insertInto(CIVICLASTMODIFIEDPROFILEIMAGE,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X32,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X256,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X14,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X64,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X128,
                        CIVICLASTMODIFIEDPROFILEIMAGE.CIVICLASTMODIFIEDORGANIZATIONID)
                        .values(userModified.organization().profileImage().x32(),
                                userModified.organization().profileImage().x256(),
                                userModified.organization().profileImage().x14(),
                                userModified.organization().profileImage().x64(),
                                userModified.organization().profileImage().x128(),
                                idLastModiefiedOrganization)
                        .execute();
            }

        }

        if (civic.lifecycleActions().lastReviewed() != null) {
            int idLastReviewed =
                    context.insertInto(CIVICLASTREVIEWED, CIVICLASTREVIEWED.TIMESTAMP, CIVICLASTREVIEWED.CIVICLIFECYCLEACTIONSID)
                            .values(civic.lifecycleActions().lastCommentedOn().timestamp(), idLifeActions)
                            .returning(CIVICLASTREVIEWED.ID)
                            .fetchOne()
                            .getValue(CIVICLASTREVIEWED.ID);

            CivicUser userLastReviewed = civic.lifecycleActions().lastReviewed().user();
            int idLastReviewedUser = context.insertInto(CIVICLASTREVIEWEDUSER,
                    CIVICLASTREVIEWEDUSER.USERNAME,
                    CIVICLASTREVIEWEDUSER.AREAOFEXPERTISE,
                    CIVICLASTREVIEWEDUSER.TWITTERHANDLE,
                    CIVICLASTREVIEWEDUSER.NAME,
                    CIVICLASTREVIEWEDUSER.BIO,
                    CIVICLASTREVIEWEDUSER.URL,
                    CIVICLASTREVIEWEDUSER.CREATEDAT,
                    CIVICLASTREVIEWEDUSER.ACCEPTEDLICENSE,
                    CIVICLASTREVIEWEDUSER.AFFILIATION,
                    CIVICLASTREVIEWEDUSER.AVATARURL,
                    CIVICLASTREVIEWEDUSER.ROLE,
                    CIVICLASTREVIEWEDUSER.FACEBOOKPROFILE,
                    CIVICLASTREVIEWEDUSER.LINKEDINPROFILE,
                    CIVICLASTREVIEWEDUSER.ORCID,
                    CIVICLASTREVIEWEDUSER.DISPLAYNAME,
                    CIVICLASTREVIEWEDUSER.LASTSEENAT,
                    CIVICLASTREVIEWEDUSER.FEATUREDEXPERT,
                    CIVICLASTREVIEWEDUSER.IDUSER,
                    CIVICLASTREVIEWEDUSER.SIGNUPCOMPLETE,
                    CIVICLASTREVIEWEDUSER.CIVICLASTREVIEWEDID)
                    .values(userLastReviewed.username(),
                            userLastReviewed.areaOfExpertise(),
                            userLastReviewed.twitterHandle(),
                            userLastReviewed.name(),
                            userLastReviewed.bio(),
                            userLastReviewed.url(),
                            userLastReviewed.createdAt(),
                            userLastReviewed.acceptedLicense(),
                            userLastReviewed.affiliation(),
                            userLastReviewed.avatarUrl(),
                            userLastReviewed.role(),
                            userLastReviewed.facebookProfile(),
                            userLastReviewed.linkedinProfile(),
                            userLastReviewed.orcid(),
                            userLastReviewed.displayName(),
                            userLastReviewed.lastSeenAt(),
                            userLastReviewed.featuredExpert(),
                            userLastReviewed.id(),
                            userLastReviewed.signupComplete(),
                            idLastReviewed)
                    .returning(CIVICLASTREVIEWEDUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTREVIEWEDUSER.ID);

            context.insertInto(CIVICLASTREVIEWEDAVATARS,
                    CIVICLASTREVIEWEDAVATARS.X32,
                    CIVICLASTREVIEWEDAVATARS.X14,
                    CIVICLASTREVIEWEDAVATARS.X64,
                    CIVICLASTREVIEWEDAVATARS.X128,
                    CIVICLASTREVIEWEDAVATARS.CIVICLASTREVIEWEDID)
                    .values(userLastReviewed.avatars().x32(),
                            userLastReviewed.avatars().x14(),
                            userLastReviewed.avatars().x64(),
                            userLastReviewed.avatars().x128(),
                            idLastReviewedUser)
                    .execute();

            int idLastReviewedOrganization = context.insertInto(CIVICLASTREVIEWEDORGANIZATION,
                    CIVICLASTREVIEWEDORGANIZATION.URL,
                    CIVICLASTREVIEWEDORGANIZATION.IDORGANIZATION,
                    CIVICLASTREVIEWEDORGANIZATION.DESCRIPTION,
                    CIVICLASTREVIEWEDORGANIZATION.NAME,
                    CIVICLASTREVIEWEDORGANIZATION.CIVICLASTREVIEWEDUSERID)
                    .values(userLastReviewed.organization().url(),
                            userLastReviewed.organization().id(),
                            userLastReviewed.organization().description(),
                            userLastReviewed.organization().name(),
                            idLastReviewedUser)
                    .returning(CIVICLASTREVIEWEDORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTREVIEWEDORGANIZATION.ID);

            context.insertInto(CIVICLASTREVIEWEDPROFILEIMAGE,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X32,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X256,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X14,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X64,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X128,
                    CIVICLASTREVIEWEDPROFILEIMAGE.CIVICLASTREVIEWEDORGANIZATIONID)
                    .values(userLastReviewed.organization().profileImage().x32(),
                            userLastReviewed.organization().profileImage().x256(),
                            userLastReviewed.organization().profileImage().x14(),
                            userLastReviewed.organization().profileImage().x64(),
                            userLastReviewed.organization().profileImage().x128(),
                            idLastReviewedOrganization)
                    .execute();

        }
    }

    private void writeKbSpecificObject(int viccEntryId, @NotNull KbSpecificObject object) {
        if (object instanceof Sage) {
            importSageinSQL(viccEntryId, object);
        }
        if (object instanceof BRCA) {
            importBRCAinSQL(viccEntryId, object);
        }
        if (object instanceof Cgi) {
            importCGIinSQL(viccEntryId, object);
        }
        if (object instanceof Jax) {
            importJAXinSQL(viccEntryId, object);
        }
        if (object instanceof JaxTrials) {
            importJaxTrialsinSQL(viccEntryId, object);
        }
        if (object instanceof Pmkb) {
            importPmkbinSQL(viccEntryId, object);
        }
        if (object instanceof Oncokb) {
            importOncokbinSQL(viccEntryId, object);
        }
        if (object instanceof MolecularMatchTrials) {
            importMolecularMatchTrials(viccEntryId, object);
        }
        if (object instanceof Civic) {
            importCivic(viccEntryId, object);
        }
        if (object instanceof MolecularMatch) {
            importMolecularMatch(viccEntryId, object);
        }
    }

    private void importMolecularMatch(int viccEntryId, @NotNull KbSpecificObject object) {
        MolecularMatch molecularMatch = (MolecularMatch) object;

        int id = context.insertInto(MOLECULARMATCH,
                MOLECULARMATCH.SCORE,
                MOLECULARMATCH.AUTOGENERATENARRATIVE,
                MOLECULARMATCH.CLINICALSIGNIFICANCE,
                MOLECULARMATCH.IDMOLECULARMATCH,
                MOLECULARMATCH.UNIQUEKEY,
                MOLECULARMATCH.CIVICVALUE,
                MOLECULARMATCH.HASHKEY,
                MOLECULARMATCH.REGULATORYBODYAPPROVED,
                MOLECULARMATCH.VERSION,
                MOLECULARMATCH.GUIDELINEBODY,
                MOLECULARMATCH.REGULATORYBODY,
                MOLECULARMATCH.CUSTOMER,
                MOLECULARMATCH.DIRECTION,
                MOLECULARMATCH.AMPCAP,
                MOLECULARMATCH.GUIDELINEVERSION,
                MOLECULARMATCH.TIER,
                MOLECULARMATCH.MVLD,
                MOLECULARMATCH.SIXTIER,
                MOLECULARMATCH.NOTHERAPYAVAILABLE,
                MOLECULARMATCH.NARRATIVE,
                MOLECULARMATCH.EXPRESSION,
                MOLECULARMATCH.BIOMARKERCLASS,
                MOLECULARMATCH.VICCENTRYID)
                .values(molecularMatch.score(),
                        molecularMatch.autoGenerateNarrative(),
                        molecularMatch.clinicalSignificance(),
                        molecularMatch.id(),
                        molecularMatch.uniqueKey(),
                        molecularMatch.civicValue(),
                        molecularMatch.hashKey(),
                        molecularMatch.regulatoryBodyApproved(),
                        molecularMatch.version(),
                        molecularMatch.guidelineBody(),
                        molecularMatch.regulatoryBody(),
                        molecularMatch.customer(),
                        molecularMatch.direction(),
                        molecularMatch.ampcap(),
                        molecularMatch.guidelineVersion(),
                        molecularMatch.tier(),
                        molecularMatch.mvld(),
                        molecularMatch.sixtier(),
                        molecularMatch.noTherapyAvailable(),
                        molecularMatch.narrative(),
                        molecularMatch.expression(),
                        molecularMatch.biomarkerClass(),
                        viccEntryId)
                .returning(MOLECULARMATCH.ID)
                .fetchOne()
                .getValue(MOLECULARMATCH.ID);

        if (molecularMatch.institution() != null) {
            for (String institution : molecularMatch.institution()) {
                context.insertInto(MOLECULARMATCHINSTUTITION,
                        MOLECULARMATCHINSTUTITION.INSTITUTION,
                        MOLECULARMATCHINSTUTITION.MOLECULARMATCHID).values(institution, id).execute();
            }
        }

        if (molecularMatch.includeGene0() != null) {
            for (String includeGene0 : molecularMatch.includeGene0()) {
                context.insertInto(MOLECULARMATCHINCLUDEGENE0,
                        MOLECULARMATCHINCLUDEGENE0.INCLUDEGENE0,
                        MOLECULARMATCHINCLUDEGENE0.MOLECULARMATCHID).values(includeGene0, id).execute();
            }
        }

        if (molecularMatch.external_id() != null) {
            for (String externalId : molecularMatch.external_id()) {
                context.insertInto(MOLECULARMATCHEXTERNALID,
                        MOLECULARMATCHEXTERNALID.EXTERNAL_ID,
                        MOLECULARMATCHEXTERNALID.MOLECULARMATCHID).values(externalId, id).execute();
            }
        }

        if (molecularMatch.includeStage0() != null) {
            for (String includeStage0 : molecularMatch.includeStage0()) {
                context.insertInto(MOLECULARMATCHINCLUDESTAGE0,
                        MOLECULARMATCHINCLUDESTAGE0.INCLUDESTAGE0,
                        MOLECULARMATCHINCLUDESTAGE0.MOLECULARMATCHID).values(includeStage0, id).execute();
            }
        }

        if (molecularMatch.includeDrug1() != null) {
            for (String includeDrug1 : molecularMatch.includeDrug1()) {
                context.insertInto(MOLECULARMATCHINCLUDEDRUG1,
                        MOLECULARMATCHINCLUDEDRUG1.INCLUDEDRUG1,
                        MOLECULARMATCHINCLUDEDRUG1.MOLECULARMATCHID).values(includeDrug1, id).execute();
            }
        }

        for (String includeCondition1 : molecularMatch.includeCondition1()) {
            context.insertInto(MOLECULARMATCHINCLUDECONDITION1,
                    MOLECULARMATCHINCLUDECONDITION1.INCLUDECONDITION1,
                    MOLECULARMATCHINCLUDECONDITION1.MOLECULARMATCHID).values(includeCondition1, id).execute();
        }

        if (molecularMatch.includeMutation1() != null) {
            for (String includeMutation1 : molecularMatch.includeMutation1()) {
                context.insertInto(MOLECULARMATCHINCLUDEMUTATION1,
                        MOLECULARMATCHINCLUDEMUTATION1.INCLUDEMUTATION1,
                        MOLECULARMATCHINCLUDEMUTATION1.MOLECULARMATCHID).values(includeMutation1, id).execute();
            }
        }

        for (String includeCondition0 : molecularMatch.includeCondition0()) {
            context.insertInto(MOLECULARMATCHINCLUDECONDITION0,
                    MOLECULARMATCHINCLUDECONDITION0.INCLUDECONDITION0,
                    MOLECULARMATCHINCLUDECONDITION0.MOLECULARMATCHID).values(includeCondition0, id).execute();
        }

        if (molecularMatch.includeMutation0() != null) {
            for (String includeMutation0 : molecularMatch.includeMutation0()) {
                context.insertInto(MOLECULARMATCHINCLUDEMUTATION0,
                        MOLECULARMATCHINCLUDEMUTATION0.INCLUDEMUTATION0,
                        MOLECULARMATCHINCLUDEMUTATION0.MOLECULARMATCHID).values(includeMutation0, id).execute();
            }
        }

        for (String criteriaMet : molecularMatch.criteriaMet()) {
            context.insertInto(MOLECULARMATCHCRITERIAMET, MOLECULARMATCHCRITERIAMET.CRITERIAMET, MOLECULARMATCHCRITERIAMET.MOLECULARMATCHID)
                    .values(criteriaMet, id)
                    .execute();
        }

        for (MolecularMatchSource source : molecularMatch.sources()) {
            context.insertInto(MOLECULARMATCHSOURCE,
                    MOLECULARMATCHSOURCE.NAME,
                    MOLECULARMATCHSOURCE.SUPPRESS,
                    MOLECULARMATCHSOURCE.PUBID,
                    MOLECULARMATCHSOURCE.SUBTYPE,
                    MOLECULARMATCHSOURCE.VALID,
                    MOLECULARMATCHSOURCE.LINK,
                    MOLECULARMATCHSOURCE.YEAR,
                    MOLECULARMATCHSOURCE.TRIALID,
                    MOLECULARMATCHSOURCE.TYPE,
                    MOLECULARMATCHSOURCE.IDSOURCE,
                    MOLECULARMATCHSOURCE.INSTITUTION,
                    MOLECULARMATCHSOURCE.TRIALPHASE,
                    MOLECULARMATCHSOURCE.FUNCTIONALCONSEQUENCE,
                    MOLECULARMATCHSOURCE.TRUSTRATING,
                    MOLECULARMATCHSOURCE.MOLECULARMATCHID)
                    .values(source.name(),
                            source.suppress(),
                            source.pubId(),
                            source.subType(),
                            source.valid(),
                            source.link(),
                            source.year(),
                            source.trialId(),
                            source.type(),
                            source.id(),
                            source.institution(),
                            source.trialPhase(),
                            source.functionalConsequence(),
                            source.trustRating(),
                            id)
                    .execute();
        }

        for (MolecularMatchTierExplanation tierExplanation : molecularMatch.tierExplanation()) {
            context.insertInto(MOLECULARMATCHTIEREXPLANATION,
                    MOLECULARMATCHTIEREXPLANATION.TIER,
                    MOLECULARMATCHTIEREXPLANATION.STEP,
                    MOLECULARMATCHTIEREXPLANATION.MESSAGE,
                    MOLECULARMATCHTIEREXPLANATION.SUCCESS,
                    MOLECULARMATCHTIEREXPLANATION.MOLECULARMATCHID)
                    .values(tierExplanation.tier(), tierExplanation.step(), tierExplanation.message(), tierExplanation.success(), id)
                    .execute();
        }

        for (MolecularMatchTherapeuticContext therapeuticContext : molecularMatch.therapeuticContext()) {
            context.insertInto(MOLECULARMATCHTHERAPEUTICCONTEXT,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.FACET,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.NAME,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.SUPPRESS,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.VALID,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.MOLECULARMATCHID)
                    .values(therapeuticContext.facet(),
                            therapeuticContext.name(),
                            therapeuticContext.suppress(),
                            therapeuticContext.valid(),
                            id)
                    .execute();
        }

        for (MolecularMatchTags tags : molecularMatch.tags()) {
            context.insertInto(MOLECULARMATCHTAGS,
                    MOLECULARMATCHTAGS.PRIORITY,
                    MOLECULARMATCHTAGS.COMPOSITEKEY,
                    MOLECULARMATCHTAGS.SUPPRESS,
                    MOLECULARMATCHTAGS.FILTERTYPE,
                    MOLECULARMATCHTAGS.TERM,
                    MOLECULARMATCHTAGS.PRIMARYVALUE,
                    MOLECULARMATCHTAGS.FACET,
                    MOLECULARMATCHTAGS.VALID,
                    MOLECULARMATCHTAGS.CUSTOM,
                    MOLECULARMATCHTAGS.ISNEW,
                    MOLECULARMATCHTAGS.GENERATEDBY,
                    MOLECULARMATCHTAGS.MANUALSUPPRESS,
                    MOLECULARMATCHTAGS.GENERATEDBYTERM,
                    MOLECULARMATCHTAGS.TRANSCRIPT,
                    MOLECULARMATCHTAGS.MOLECULARMATCHID)
                    .values(tags.priority(),
                            tags.compositeKey(),
                            tags.suppress(),
                            tags.filterType(),
                            tags.term(),
                            tags.primary(),
                            tags.facet(),
                            tags.valid(),
                            tags.custom(),
                            tags.isNew(),
                            tags.generatedBy(),
                            tags.manualSuppress(),
                            tags.generatedByTerm(),
                            tags.transcript(),
                            id)
                    .execute();
        }

        for (MolecularMatchVariantInfo variantInfo : molecularMatch.variantInfo()) {
            int idVariantInfo = context.insertInto(MOLECULARMATCHVARIANTINFO,
                    MOLECULARMATCHVARIANTINFO.CLASSIFICATION,
                    MOLECULARMATCHVARIANTINFO.NAME,
                    MOLECULARMATCHVARIANTINFO.GENEFUSIONPARTNER,
                    MOLECULARMATCHVARIANTINFO.COSMIC_ID,
                    MOLECULARMATCHVARIANTINFO.GENE,
                    MOLECULARMATCHVARIANTINFO.TRANSCRIPT,
                    MOLECULARMATCHVARIANTINFO.POPFREQMAX,
                    MOLECULARMATCHVARIANTINFO.MOLECULARMATCHID)
                    .values(variantInfo.classification(),
                            variantInfo.name(),
                            variantInfo.geneFusionPartner(),
                            variantInfo.COSMIC_ID(),
                            variantInfo.gene(),
                            variantInfo.transcript(),
                            variantInfo.popFreqMax(),
                            id)
                    .returning(MOLECULARMATCHVARIANTINFO.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHVARIANTINFO.ID);

            for (String consequences : variantInfo.consequences()) {
                context.insertInto(MOLECULARMATCHVARIANTINFOCONSEQUENCES,
                        MOLECULARMATCHVARIANTINFOCONSEQUENCES.CONSEQUENCES,
                        MOLECULARMATCHVARIANTINFOCONSEQUENCES.MOLECULARMATCHVARIANTINFOID).values(consequences, idVariantInfo).execute();
            }

            for (MolecularMatchFusions fusions : variantInfo.fusions()) {
                context.insertInto(MOLECULARMATCHVARIANTINFOFUSIONS,
                        MOLECULARMATCHVARIANTINFOFUSIONS.REFERENCEGENOME,
                        MOLECULARMATCHVARIANTINFOFUSIONS.LBPWREP,
                        MOLECULARMATCHVARIANTINFOFUSIONS.RBPWREP,
                        MOLECULARMATCHVARIANTINFOFUSIONS.EXONNUMBER,
                        MOLECULARMATCHVARIANTINFOFUSIONS.CHR,
                        MOLECULARMATCHVARIANTINFOFUSIONS.RBPWLEP,
                        MOLECULARMATCHVARIANTINFOFUSIONS.INTRONNUMBER,
                        MOLECULARMATCHVARIANTINFOFUSIONS.LBPWLEP,
                        MOLECULARMATCHVARIANTINFOFUSIONS.MOLECULARMATCHVARIANTINFOID)
                        .values(fusions.referenceGenome(),
                                fusions.LBPWREP(),
                                fusions.RBPWREP(),
                                fusions.exonNumber(),
                                fusions.chr(),
                                fusions.RBPWLEP(),
                                fusions.intronNumber(),
                                fusions.RBPWLEP(),
                                idVariantInfo)
                        .execute();
            }

            for (MolecularMatchLocations locations : variantInfo.locations()) {
                int Idlocation = context.insertInto(MOLECULARMATCHVARIANTINFOLOCATIONS,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.AMINOACIDCHANGE,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.INTRONNUMBER,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.STOP,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.START,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.CHR,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.STRAND,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.ALT,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.REFERENCEGENOME,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.REF,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.CDNA,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.MOLECULARMATCHVARIANTINFOID)
                        .values(locations.aminoAcidChange(),
                                locations.intronNumber(),
                                locations.stop(),
                                locations.start(),
                                locations.chr(),
                                locations.strand(),
                                locations.alt(),
                                locations.referenceGenome(),
                                locations.ref(),
                                locations.cdna(),
                                idVariantInfo)
                        .returning(MOLECULARMATCHVARIANTINFOLOCATIONS.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHVARIANTINFOLOCATIONS.ID);

                if (locations.exonNumber() != null) {
                    for (String exonNumber : locations.exonNumber()) {
                        context.insertInto(MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER,
                                MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER.EXONNUMBER,
                                MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER.MOLECULARMATCHVARIANTINFOLOCATIONSID)
                                .values(exonNumber, Idlocation)
                                .execute();

                    }
                }
            }
        }

        for (MolecularMatchClassification classification : molecularMatch.classification()) {
            int idClassification = context.insertInto(MOLECULARMATCHCLASSIFICATION,
                    MOLECULARMATCHCLASSIFICATION.CLASSIFICATION,
                    MOLECULARMATCHCLASSIFICATION.CLASSIFICATIONOVERRIDE,
                    MOLECULARMATCHCLASSIFICATION.GENESYMBOL,
                    MOLECULARMATCHCLASSIFICATION.DESCRIPTION,
                    MOLECULARMATCHCLASSIFICATION.PRIORITY,
                    MOLECULARMATCHCLASSIFICATION.EXPANDGENESEARCH,
                    MOLECULARMATCHCLASSIFICATION.DRUGSEXPERIMENTALCOUNT,
                    MOLECULARMATCHCLASSIFICATION.DRUGSAPPROVEDOFFLABELCOUNT,
                    MOLECULARMATCHCLASSIFICATION.COPYNUMBERTYPE,
                    MOLECULARMATCHCLASSIFICATION.PUBLICATIONCOUNT,
                    MOLECULARMATCHCLASSIFICATION.TRANSCRIPT,
                    MOLECULARMATCHCLASSIFICATION.NAME,
                    MOLECULARMATCHCLASSIFICATION.ROOTTERM,
                    MOLECULARMATCHCLASSIFICATION.DRUGSAPPROVEDONLABELCOUNT,
                    MOLECULARMATCHCLASSIFICATION.TRIALCOUNT,
                    MOLECULARMATCHCLASSIFICATION.ALIAS,
                    MOLECULARMATCHCLASSIFICATION.MOLECULARMATCHID)
                    .values(classification.classification(),
                            classification.classificationOverride(),
                            classification.geneSymbol(),
                            classification.description(),
                            classification.priority(),
                            classification.expandGeneSearch(),
                            classification.drugsExperimentalCount(),
                            classification.drugsApprovedOffLabelCount(),
                            classification.copyNumberType(),
                            classification.publicationCount(),
                            classification.transcript(),
                            classification.name(),
                            classification.rootTerm(),
                            classification.drugsApprovedOnLabelCount(),
                            classification.trialCount(),
                            classification.alias(),
                            id)
                    .returning(MOLECULARMATCHCLASSIFICATION.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHCLASSIFICATION.ID);

            if (classification.end() != null) {
                for (String end : classification.end()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONEND,
                            MOLECULARMATCHCLASSIFICATIONEND.END,
                            MOLECULARMATCHCLASSIFICATIONEND.MOLECULARMATCHCLASSIFICATIONID).values(end, idClassification).execute();
                }
            }

            if (classification.start() != null) {
                for (String start : classification.start()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONSTART,
                            MOLECULARMATCHCLASSIFICATIONSTART.START,
                            MOLECULARMATCHCLASSIFICATIONSTART.MOLECULARMATCHCLASSIFICATIONID).values(start, idClassification).execute();
                }
            }

            if (classification.chr() != null) {
                for (String chr : classification.chr()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONCHR,
                            MOLECULARMATCHCLASSIFICATIONCHR.CHROMOSOME,
                            MOLECULARMATCHCLASSIFICATIONCHR.MOLECULARMATCHCLASSIFICATIONID).values(chr, idClassification).execute();
                }
            }

            if (classification.pathology() != null) {
                for (String pathology : classification.pathology()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONPATHOLOGY,
                            MOLECULARMATCHCLASSIFICATIONPATHOLOGY.PATHOLOGY,
                            MOLECULARMATCHCLASSIFICATIONPATHOLOGY.MOLECULARMATCHCLASSIFICATIONID)
                            .values(pathology, idClassification)
                            .execute();
                }
            }

            if (classification.ref() != null) {
                for (String ref : classification.ref()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONREF,
                            MOLECULARMATCHCLASSIFICATIONREF.REF,
                            MOLECULARMATCHCLASSIFICATIONREF.MOLECULARMATCHCLASSIFICATIONID).values(ref, idClassification).execute();
                }
            }

            if (classification.exon() != null) {
                for (String exon : classification.exon()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONEXON,
                            MOLECULARMATCHCLASSIFICATIONEXON.EXON,
                            MOLECULARMATCHCLASSIFICATIONEXON.MOLECULARMATCHCLASSIFICATIONID).values(exon, idClassification).execute();
                }
            }

            if (classification.alt() != null) {
                for (String alt : classification.alt()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONALT,
                            MOLECULARMATCHCLASSIFICATIONALT.ALT,
                            MOLECULARMATCHCLASSIFICATIONALT.MOLECULARMATCHCLASSIFICATIONID).values(alt, idClassification).execute();
                }
            }

            if (classification.sources() != null) {
                for (String source : classification.sources()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONSOURCES,
                            MOLECULARMATCHCLASSIFICATIONSOURCES.SOURCE,
                            MOLECULARMATCHCLASSIFICATIONSOURCES.MOLECULARMATCHCLASSIFICATIONID).values(source, idClassification).execute();
                }
            }

            if (classification.transcripts() != null) {
                for (String transcript : classification.transcripts()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS,
                            MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS.TRANSCRIPT,
                            MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS.MOLECULARMATCHCLASSIFICATIONID)
                            .values(transcript, idClassification)
                            .execute();
                }
            }

            if (classification.COSMIC_ID() != null) {
                for (String cosmicId : classification.COSMIC_ID()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONCOSMICID,
                            MOLECULARMATCHCLASSIFICATIONCOSMICID.COSMIC_ID,
                            MOLECULARMATCHCLASSIFICATIONCOSMICID.MOLECULARMATCHCLASSIFICATIONID)
                            .values(cosmicId, idClassification)
                            .execute();
                }
            }

            if (classification.dbSNP() != null) {
                for (String dbSNP : classification.dbSNP()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONDBSNP,
                            MOLECULARMATCHCLASSIFICATIONDBSNP.DBSNP,
                            MOLECULARMATCHCLASSIFICATIONDBSNP.MOLECULARMATCHCLASSIFICATIONID).values(dbSNP, idClassification).execute();
                }
            }

            if (classification.popFreqMax() != null) {
                for (String popFreqMax : classification.popFreqMax()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONPOPFREQMAX,
                            MOLECULARMATCHCLASSIFICATIONPOPFREQMAX.POPFREQMAX,
                            MOLECULARMATCHCLASSIFICATIONPOPFREQMAX.MOLECULARMATCHCLASSIFICATIONID)
                            .values(popFreqMax, idClassification)
                            .execute();
                }
            }

            if (classification.exonicFunc() != null) {
                for (String exonicFunc : classification.exonicFunc()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONEXONICFUNC,
                            MOLECULARMATCHCLASSIFICATIONEXONICFUNC.EXONICFUNC,
                            MOLECULARMATCHCLASSIFICATIONEXONICFUNC.MOLECULARMATCHCLASSIFICATIONID)
                            .values(exonicFunc, idClassification)
                            .execute();
                }
            }

            if (classification.NucleotideChange() != null) {
                for (String nucleotideChange : classification.NucleotideChange()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE,
                            MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE.NUCLEOTIDECHANGE,
                            MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE.MOLECULARMATCHCLASSIFICATIONID)
                            .values(nucleotideChange, idClassification)
                            .execute();
                }
            }

            if (classification.parents() != null) {
                for (MolecularMatchParents parents : classification.parents()) {
                    int idParents = context.insertInto(MOLECULARMATCHCLASSIFICATIONPARENTS,
                            MOLECULARMATCHCLASSIFICATIONPARENTS.TYPE,
                            MOLECULARMATCHCLASSIFICATIONPARENTS.NAME,
                            MOLECULARMATCHCLASSIFICATIONPARENTS.ACTIONABLEPARENT,
                            MOLECULARMATCHCLASSIFICATIONPARENTS.MOLECULARMATCHCLASSIFICATIONID)
                            .values(parents.type(), parents.name(), parents.actionableParent(), idClassification)
                            .returning(MOLECULARMATCHCLASSIFICATIONPARENTS.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHCLASSIFICATIONPARENTS.ID);

                    for (String transcripts : parents.transcripts()) {
                        context.insertInto(MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS,
                                MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS.TRANSCRIPT,
                                MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS.MOLECULARMATCHCLASSIFICATIONPARENTSID)
                                .values(transcripts, idParents)
                                .execute();
                    }
                }
            }
        }

        for (MolecularMatchCriteriaUnmet criteriaUnmet : molecularMatch.criteriaUnmet()) {
            context.insertInto(MOLECULARMATCHCRITERIAUNMET,
                    MOLECULARMATCHCRITERIAUNMET.PRIORITY,
                    MOLECULARMATCHCRITERIAUNMET.COMPOSITEKEY,
                    MOLECULARMATCHCRITERIAUNMET.ISNEW,
                    MOLECULARMATCHCRITERIAUNMET.GENERATEDBY,
                    MOLECULARMATCHCRITERIAUNMET.MANUALSUPPRESS,
                    MOLECULARMATCHCRITERIAUNMET.GENERATEDBYTERM,
                    MOLECULARMATCHCRITERIAUNMET.SUPPRESS,
                    MOLECULARMATCHCRITERIAUNMET.FILTERTYPE,
                    MOLECULARMATCHCRITERIAUNMET.TERM,
                    MOLECULARMATCHCRITERIAUNMET.PRIMARYVALUE,
                    MOLECULARMATCHCRITERIAUNMET.FACET,
                    MOLECULARMATCHCRITERIAUNMET.VALID,
                    MOLECULARMATCHCRITERIAUNMET.CUSTOM,
                    MOLECULARMATCHCRITERIAUNMET.TRANSCRIPT,
                    MOLECULARMATCHCRITERIAUNMET.MOLECULARMATCHID)
                    .values(criteriaUnmet.priority(),
                            criteriaUnmet.compositeKey(),
                            criteriaUnmet.isNew(),
                            criteriaUnmet.generatedBy(),
                            criteriaUnmet.manualSuppress(),
                            criteriaUnmet.generatedByTerm(),
                            criteriaUnmet.suppress(),
                            criteriaUnmet.filterType(),
                            criteriaUnmet.term(),
                            criteriaUnmet.primary(),
                            criteriaUnmet.facet(),
                            criteriaUnmet.valid(),
                            criteriaUnmet.custom(),
                            criteriaUnmet.transcript(),
                            id)
                    .execute();
        }

        for (MolecularMatchPrevalence prevelance : molecularMatch.prevalence()) {
            context.insertInto(MOLECULARMATCHPREFELANCE,
                    MOLECULARMATCHPREFELANCE.COUNT,
                    MOLECULARMATCHPREFELANCE.PERCENT,
                    MOLECULARMATCHPREFELANCE.STUDYID,
                    MOLECULARMATCHPREFELANCE.SAMPLES,
                    MOLECULARMATCHPREFELANCE.MOLECULAR,
                    MOLECULARMATCHPREFELANCE.CONDITIONVALUE,
                    MOLECULARMATCHPREFELANCE.MOLECULARMATCHID)
                    .values(prevelance.count(),
                            prevelance.percent(),
                            prevelance.studyId(),
                            prevelance.samples(),
                            prevelance.molecular(),
                            prevelance.condition(),
                            id)
                    .execute();
        }

        for (MolecularMatchMutations mutations : molecularMatch.mutations()) {
            int idMutations = context.insertInto(MOLECULARMATCHMUTATIONS,
                    MOLECULARMATCHMUTATIONS.LONGESTTRANSCRIPT,
                    MOLECULARMATCHMUTATIONS.DESCRIPTION,
                    MOLECULARMATCHMUTATIONS.SRC,
                    MOLECULARMATCHMUTATIONS.UNIPROTTRANSCRIPT,
                    MOLECULARMATCHMUTATIONS.TRANSCRIPTRECOGNIZED,
                    MOLECULARMATCHMUTATIONS.GENESYMBOL,
                    MOLECULARMATCHMUTATIONS.TRANSCRIPT,
                    MOLECULARMATCHMUTATIONS.IDMUTATIONS,
                    MOLECULARMATCHMUTATIONS.NAME,
                    MOLECULARMATCHMUTATIONS.MOLECULARMATCHID)
                    .values(mutations.longestTranscript(),
                            mutations.description(),
                            mutations.src(),
                            mutations.uniprotTranscript(),
                            mutations.transcriptRecognized(),
                            mutations.geneSymbol(),
                            mutations.transcript(),
                            mutations.id(),
                            mutations.name(),
                            id)
                    .returning(MOLECULARMATCHMUTATIONS.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHMUTATIONS.ID);

            for (String mutationType : mutations.mutationType()) {
                context.insertInto(MOLECULARMATCHMUTATIONSMUTATIONTYPE,
                        MOLECULARMATCHMUTATIONSMUTATIONTYPE.MUTATIONTYPE,
                        MOLECULARMATCHMUTATIONSMUTATIONTYPE.MOLECULARMATCHMUTATIONSID).values(mutationType, idMutations).execute();
            }

            for (String source : mutations.sources()) {
                context.insertInto(MOLECULARMATCHMUTATIONSSOURCE,
                        MOLECULARMATCHMUTATIONSSOURCE.SOURCE,
                        MOLECULARMATCHMUTATIONSSOURCE.MOLECULARMATCHMUTATIONSID).values(source, idMutations).execute();
            }

            for (String synonyms : mutations.synonyms()) {
                context.insertInto(MOLECULARMATCHMUTATIONSSYNONYMS,
                        MOLECULARMATCHMUTATIONSSYNONYMS.SYNONYMS,
                        MOLECULARMATCHMUTATIONSSYNONYMS.MOLECULARMATCHMUTATIONSID).values(synonyms, idMutations).execute();
            }

            for (String pathology : mutations.pathology()) {
                context.insertInto(MOLECULARMATCHMUTATIONSPATHOLOGY,
                        MOLECULARMATCHMUTATIONSPATHOLOGY.PATHOLOGY,
                        MOLECULARMATCHMUTATIONSPATHOLOGY.MOLECULARMATCHMUTATIONSID).values(pathology, idMutations).execute();
            }

            for (String cDNA : mutations.cDNA()) {
                context.insertInto(MOLECULARMATCHMUTATIONSCDNA,
                        MOLECULARMATCHMUTATIONSCDNA.CDNA,
                        MOLECULARMATCHMUTATIONSCDNA.MOLECULARMATCHMUTATIONSID).values(cDNA, idMutations).execute();
            }

            for (MolecularMatchTranscriptConsequence transcriptConsequence : mutations.transcriptConsequence()) {
                int idTranscriptConsequence = context.insertInto(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.AMINOACIDCHANGE,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.COMPOSITEKEY,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.INTRONNUMBER,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.SUPPRESS,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.STOP,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.CUSTOM,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.START,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.CHR,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.STRAND,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.VALIDATED,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.TRANSCRIPT,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.CDNA,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.REFERENCEGENOME,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.REF,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.ALT,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.MOLECULARMATCHMUTATIONSID)
                        .values(transcriptConsequence.aminoAcidChange(),
                                transcriptConsequence.compositeKey(),
                                transcriptConsequence.intronNumber(),
                                transcriptConsequence.suppress(),
                                transcriptConsequence.stop(),
                                transcriptConsequence.custom(),
                                transcriptConsequence.start(),
                                transcriptConsequence.chr(),
                                transcriptConsequence.strand(),
                                transcriptConsequence.validated(),
                                transcriptConsequence.transcript(),
                                transcriptConsequence.cdna(),
                                transcriptConsequence.referenceGenome(),
                                transcriptConsequence.ref(),
                                transcriptConsequence.alt(),
                                idMutations)
                        .returning(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.ID);

                if (transcriptConsequence.exonNumber() != null) {
                    for (String exonNumber : transcriptConsequence.exonNumber()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER,
                                MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER.EXONNUMBER,
                                MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER.MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESID)
                                .values(exonNumber, idTranscriptConsequence)
                                .execute();
                    }
                }
            }

            for (MolecularMatchParents parents : mutations.parents()) {
                int idParents = context.insertInto(MOLECULARMATCHMUTATIONSPARENTS,
                        MOLECULARMATCHMUTATIONSPARENTS.TYPE,
                        MOLECULARMATCHMUTATIONSPARENTS.NAME,
                        MOLECULARMATCHMUTATIONSPARENTS.ACTIONABLEPARENT,
                        MOLECULARMATCHMUTATIONSPARENTS.MOLECULARMATCHMUTATIONSID)
                        .values(parents.type(), parents.name(), parents.actionableParent(), idMutations)
                        .returning(MOLECULARMATCHMUTATIONSPARENTS.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSPARENTS.ID);

                for (String transcripts : parents.transcripts()) {
                    context.insertInto(MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT,
                            MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT.MOLECULARMATCHMUTATIONSPARENTSID)
                            .values(transcripts, idParents)
                            .execute();
                }
            }

            if (mutations.wgsaData() != null) {
                for (MolecularMatchWGSadataLocation wgsaDataLocation : mutations.wgsaData()) {
                    int idWGSData = context.insertInto(MOLECULARMATCHMUTATIONSWGSDATA,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXONICFUNC,
                            MOLECULARMATCHMUTATIONSWGSDATA.DBSNP,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_NFE,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_FIN,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_ALL,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_SAS,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_EAS,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_AFR,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_SAS,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_EAS,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_AMR,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_AFR,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_FREQ,
                            MOLECULARMATCHMUTATIONSWGSDATA.END,
                            MOLECULARMATCHMUTATIONSWGSDATA.START,
                            MOLECULARMATCHMUTATIONSWGSDATA.SIPHY_29WAY_LOGODDS,
                            MOLECULARMATCHMUTATIONSWGSDATA.REF,
                            MOLECULARMATCHMUTATIONSWGSDATA.GERP_RS,
                            MOLECULARMATCHMUTATIONSWGSDATA.FATHMM,
                            MOLECULARMATCHMUTATIONSWGSDATA.NUCLEOTIDECHANGE,
                            MOLECULARMATCHMUTATIONSWGSDATA.PHYLOP100WAY_VERTEBRATE,
                            MOLECULARMATCHMUTATIONSWGSDATA.FUNC,
                            MOLECULARMATCHMUTATIONSWGSDATA.GWAS_PUBMED,
                            MOLECULARMATCHMUTATIONSWGSDATA.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONSWGSDATA.ESP6500SI_AA,
                            MOLECULARMATCHMUTATIONSWGSDATA.ESP6500SI_EA,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_EUR,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_AMR,
                            MOLECULARMATCHMUTATIONSWGSDATA.CHR_START_REF_ALT,
                            MOLECULARMATCHMUTATIONSWGSDATA.AA,
                            MOLECULARMATCHMUTATIONSWGSDATA.POPFREQMAX,
                            MOLECULARMATCHMUTATIONSWGSDATA.FATHMM_PRED,
                            MOLECULARMATCHMUTATIONSWGSDATA.WGRNA,
                            MOLECULARMATCHMUTATIONSWGSDATA.PHYLOP46WAY_PLACENTAL,
                            MOLECULARMATCHMUTATIONSWGSDATA.KEYVALUE,
                            MOLECULARMATCHMUTATIONSWGSDATA.TARGETSCANS,
                            MOLECULARMATCHMUTATIONSWGSDATA.CHR,
                            MOLECULARMATCHMUTATIONSWGSDATA.COSMIC_ID,
                            MOLECULARMATCHMUTATIONSWGSDATA.ALT,
                            MOLECULARMATCHMUTATIONSWGSDATA.GWAS_DIS,
                            MOLECULARMATCHMUTATIONSWGSDATA.GWAS_SNP,
                            MOLECULARMATCHMUTATIONSWGSDATA.MOLECULARMATCHMUTATIONSID)
                            .values(wgsaDataLocation.ExonicFunc(),
                                    wgsaDataLocation.dbSNP(),
                                    wgsaDataLocation.ExAC_NFE(),
                                    wgsaDataLocation.ExAC_FIN(),
                                    wgsaDataLocation.G1000_ALL(),
                                    wgsaDataLocation.G1000_SAS(),
                                    wgsaDataLocation.G1000_EAS(),
                                    wgsaDataLocation.G1000_AFR(),
                                    wgsaDataLocation.ExAC_SAS(),
                                    wgsaDataLocation.ExAC_EAS(),
                                    wgsaDataLocation.ExAC_AMR(),
                                    wgsaDataLocation.ExAC_AFR(),
                                    wgsaDataLocation.ExAC_Freq(),
                                    wgsaDataLocation.End(),
                                    wgsaDataLocation.Start(),
                                    wgsaDataLocation.SiPhy_29way_logOdds(),
                                    wgsaDataLocation.Ref(),
                                    wgsaDataLocation.GERP_RS(),
                                    wgsaDataLocation.FATHMM(),
                                    wgsaDataLocation.NucleotideChange(),
                                    wgsaDataLocation.phyloP100way_vertebrate(),
                                    wgsaDataLocation.Func(),
                                    wgsaDataLocation.GWAS_PUBMED(),
                                    wgsaDataLocation.Transcript(),
                                    wgsaDataLocation.ESP6500si_AA(),
                                    wgsaDataLocation.ESP6500si_EA(),
                                    wgsaDataLocation.G1000_EUR(),
                                    wgsaDataLocation.G1000_AMR(),
                                    wgsaDataLocation.Chr_Start_Ref_Alt(),
                                    wgsaDataLocation.AA(),
                                    wgsaDataLocation.PopFreqMax(),
                                    wgsaDataLocation.FATHMM_Pred(),
                                    wgsaDataLocation.wgRna(),
                                    wgsaDataLocation.phyloP46way_placental(),
                                    wgsaDataLocation.key(),
                                    wgsaDataLocation.targetScanS(),
                                    wgsaDataLocation.Chr(),
                                    wgsaDataLocation.COSMIC_ID(),
                                    wgsaDataLocation.alt(),
                                    wgsaDataLocation.GWAS_DIS(),
                                    wgsaDataLocation.GWAS_SNP(),
                                    idMutations)
                            .returning(MOLECULARMATCHMUTATIONSWGSDATA.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHMUTATIONSWGSDATA.ID);

                    if (wgsaDataLocation.ClinVar_DIS() != null) {
                        for (String clinvarDIS : wgsaDataLocation.ClinVar_DIS()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS.CLINVARDIS,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS.MOLECULARMATCHMUTATIONSWGSDATAID)
                                    .values(clinvarDIS, idWGSData)
                                    .execute();
                        }
                    }

                    if (wgsaDataLocation.ClinVar_SIG() != null) {
                        for (String clinvarSig : wgsaDataLocation.ClinVar_SIG()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG.CLINVARSIG,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG.MOLECULARMATCHMUTATIONSWGSDATAID)
                                    .values(clinvarSig, idWGSData)
                                    .execute();
                        }
                    }

                    if (wgsaDataLocation.ClinVar_STATUS() != null) {
                        for (String clinvarStatus : wgsaDataLocation.ClinVar_STATUS()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS.CLINVARSTATUS,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS.MOLECULARMATCHMUTATIONSWGSDATAID)
                                    .values(clinvarStatus, idWGSData)
                                    .execute();
                        }
                    }

                    if (wgsaDataLocation.ClinVar_DBID() != null) {
                        for (String clinvarDBID : wgsaDataLocation.ClinVar_DBID()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID.CLINVARDBID,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID.MOLECULARMATCHMUTATIONSWGSDATAID)
                                    .values(clinvarDBID, idWGSData)
                                    .execute();
                        }
                    }

                    for (String fullAA : wgsaDataLocation.FullAA()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSWGSDATAFULLAA,
                                MOLECULARMATCHMUTATIONSWGSDATAFULLAA.MOLECULARMATCHMUTATIONSWGSDATAFULLAA.FULLAA,
                                MOLECULARMATCHMUTATIONSWGSDATAFULLAA.MOLECULARMATCHMUTATIONSWGSDATAID).values(fullAA, idWGSData).execute();
                    }

                    for (String gene : wgsaDataLocation.Gene()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSWGSDATAGENE,
                                MOLECULARMATCHMUTATIONSWGSDATAGENE.GENE,
                                MOLECULARMATCHMUTATIONSWGSDATAGENE.MOLECULARMATCHMUTATIONSWGSDATAID).values(gene, idWGSData).execute();
                    }
                }
            }

            if (mutations.wgsaMap() != null) {
                for (MolecularMatchWGSaMap wgSaMap : mutations.wgsaMap()) {
                    int idMap = context.insertInto(MOLECULARMATCHMUTATIONSWGSMAP,
                            MOLECULARMATCHMUTATIONSWGSMAP.AA,
                            MOLECULARMATCHMUTATIONSWGSMAP.NAME,
                            MOLECULARMATCHMUTATIONSWGSMAP.GRCH37_CHR_START_REF_ALT,
                            MOLECULARMATCHMUTATIONSWGSMAP.NUCLEOTIDECHANGE,
                            MOLECULARMATCHMUTATIONSWGSMAP.EXON,
                            MOLECULARMATCHMUTATIONSWGSMAP.GENE,
                            MOLECULARMATCHMUTATIONSWGSMAP.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONSWGSMAP.MOLECULARMATCHMUTATIONSID)
                            .values(wgSaMap.AA(),
                                    wgSaMap.name(),
                                    wgSaMap.GRCh37_Chr_Start_Ref_Alt(),
                                    wgSaMap.NucleotideChange(),
                                    wgSaMap.Exon(),
                                    wgSaMap.Gene(),
                                    wgSaMap.Transcript(),
                                    idMutations)
                            .returning(MOLECULARMATCHMUTATIONSWGSMAP.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHMUTATIONSWGSMAP.ID);

                    for (String synonyms : wgSaMap.Synonyms()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS,
                                MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS.SYNONYMS,
                                MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS.MOLECULARMATCHMUTATIONSWGSMAPID).values(synonyms, idMap).execute();
                    }

                    for (String protCoords : wgSaMap.ProtCoords()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS,
                                MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS.PROTCOORDS,
                                MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS.MOLECULARMATCHMUTATIONSWGSMAPID)
                                .values(protCoords, idMap)
                                .execute();
                    }

                }
            }

            for (MolecularMatchGRch37Location gRch37Location : mutations.gRch37Location()) {
                int idLocation = context.insertInto(MOLECULARMATCHMUTATIONSGRCH37LOCATION,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.COMPOSITEKEY,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.REF,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.STOP,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.START,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.CHR,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.ALT,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.VALIDATED,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.STRAND,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.MOLECULARMATCHMUTATIONSID)
                        .values(gRch37Location.compositeKey(),
                                gRch37Location.ref(),
                                gRch37Location.stop(),
                                gRch37Location.start(),
                                gRch37Location.chr(),
                                gRch37Location.alt(),
                                gRch37Location.validated(),
                                gRch37Location.strand(),
                                idMutations)
                        .returning(MOLECULARMATCHMUTATIONSGRCH37LOCATION.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSGRCH37LOCATION.ID);

                for (MolecularMatchTranscriptConsequencesGRCH37 transcriptConsequencesGRCH37 : gRch37Location.transcriptConsequences()) {
                    int idConsequences = context.insertInto(MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.AMINOACIDCHANGE,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.INTRONNUMBER,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.CDNA,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.MOLECULARMATCHMUTATIONSGRCH37LOCATIONID)
                            .values(transcriptConsequencesGRCH37.aminoAcidChange(),
                                    transcriptConsequencesGRCH37.intronNumber(),
                                    transcriptConsequencesGRCH37.transcript(),
                                    transcriptConsequencesGRCH37.cdna(),
                                    idLocation)
                            .returning(MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.ID);

                    for (String txSites : transcriptConsequencesGRCH37.txSites()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES,
                                MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES.TXSITES,
                                MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES.MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCESID)
                                .values(txSites, idConsequences)
                                .execute();
                    }

                    if (transcriptConsequencesGRCH37.exonNumber() != null) {
                        for (String exonNumber : transcriptConsequencesGRCH37.exonNumber()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER,
                                    MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER.EXONNUMBER,
                                    MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER.MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCESID)
                                    .values(exonNumber, idConsequences)
                                    .execute();
                        }
                    }
                }
            }

            if (mutations.fusionData() != null) {
                for (MolecularMatchFusionData fusions : mutations.fusionData()) {
                    int idFusions = context.insertInto(MOLECULARMATCHMUTATIONSFUSION,
                            MOLECULARMATCHMUTATIONSFUSION.SYNONYM,
                            MOLECULARMATCHMUTATIONSFUSION.SOURCE,
                            MOLECULARMATCHMUTATIONSFUSION.PAPER,
                            MOLECULARMATCHMUTATIONSFUSION.MOLECULARMATCHMUTATIONSID)
                            .values(fusions.synonym(), fusions.source(), fusions.Paper(), idMutations)
                            .returning(MOLECULARMATCHMUTATIONSFUSION.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHMUTATIONSFUSION.ID);

                    if (fusions.Bchr() != null) {
                        for (String Bchr : fusions.Bchr()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBCHR,
                                    MOLECULARMATCHMUTATIONSFUSIONBCHR.BRCHR,
                                    MOLECULARMATCHMUTATIONSFUSIONBCHR.MOLECULARMATCHMUTATIONSFUSIONID).values(Bchr, idFusions).execute();
                        }
                    }

                    if (fusions.Agene() != null) {
                        for (String agene : fusions.Agene()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONAGENE,
                                    MOLECULARMATCHMUTATIONSFUSIONAGENE.AGENE,
                                    MOLECULARMATCHMUTATIONSFUSIONAGENE.MOLECULARMATCHMUTATIONSFUSIONID).values(agene, idFusions).execute();
                        }
                    }

                    if (fusions.Btx() != null) {
                        for (String btx : fusions.Btx()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBTX,
                                    MOLECULARMATCHMUTATIONSFUSIONBTX.BTX,
                                    MOLECULARMATCHMUTATIONSFUSIONBTX.MOLECULARMATCHMUTATIONSFUSIONID).values(btx, idFusions).execute();
                        }
                    }

                    if (fusions.Achr() != null) {
                        for (String achr : fusions.Achr()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONATX,
                                    MOLECULARMATCHMUTATIONSFUSIONATX.ATX,
                                    MOLECULARMATCHMUTATIONSFUSIONATX.MOLECULARMATCHMUTATIONSFUSIONID).values(achr, idFusions).execute();
                        }
                    }

                    if (fusions.ins() != null) {
                        for (String ins : fusions.ins()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONINS,
                                    MOLECULARMATCHMUTATIONSFUSIONINS.INS,
                                    MOLECULARMATCHMUTATIONSFUSIONINS.MOLECULARMATCHMUTATIONSFUSIONID).values(ins, idFusions).execute();
                        }
                    }

                    if (fusions.Bgene() != null) {
                        for (String bgene : fusions.Bgene()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBGENE,
                                    MOLECULARMATCHMUTATIONSFUSIONBGENE.BGENE,
                                    MOLECULARMATCHMUTATIONSFUSIONBGENE.MOLECULARMATCHMUTATIONSFUSIONID).values(bgene, idFusions).execute();
                        }
                    }

                    if (fusions.Acoord() != null) {
                        for (String Acoord : fusions.Acoord()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONACOORD,
                                    MOLECULARMATCHMUTATIONSFUSIONACOORD.ACOORD,
                                    MOLECULARMATCHMUTATIONSFUSIONACOORD.MOLECULARMATCHMUTATIONSFUSIONID)
                                    .values(Acoord, idFusions)
                                    .execute();
                        }
                    }

                    if (fusions.Bori() != null) {
                        for (String bori : fusions.Bori()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBORI,
                                    MOLECULARMATCHMUTATIONSFUSIONBORI.BORI,
                                    MOLECULARMATCHMUTATIONSFUSIONBORI.MOLECULARMATCHMUTATIONSFUSIONID).values(bori, idFusions).execute();
                        }
                    }

                    if (fusions.Aband() != null) {
                        for (String aband : fusions.Aband()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONABAND,
                                    MOLECULARMATCHMUTATIONSFUSIONABAND.ABAND,
                                    MOLECULARMATCHMUTATIONSFUSIONABAND.MOLECULARMATCHMUTATIONSFUSIONID).values(aband, idFusions).execute();
                        }
                    }

                    if (fusions.Bband() != null) {
                        for (String bband : fusions.Bband()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBBAND,
                                    MOLECULARMATCHMUTATIONSFUSIONBBAND.BBAND,
                                    MOLECULARMATCHMUTATIONSFUSIONBBAND.MOLECULARMATCHMUTATIONSFUSIONID).values(bband, idFusions).execute();
                        }
                    }

                    if (fusions.Aori() != null) {
                        for (String aori : fusions.Aori()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONAORI,
                                    MOLECULARMATCHMUTATIONSFUSIONAORI.AORI,
                                    MOLECULARMATCHMUTATIONSFUSIONAORI.MOLECULARMATCHMUTATIONSFUSIONID).values(aori, idFusions).execute();
                        }
                    }

                    if (fusions.Atx() != null) {
                        for (String atx : fusions.Atx()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONATX,
                                    MOLECULARMATCHMUTATIONSFUSIONATX.ATX,
                                    MOLECULARMATCHMUTATIONSFUSIONATX.MOLECULARMATCHMUTATIONSFUSIONID).values(atx, idFusions).execute();
                        }
                    }

                    if (fusions.Bcoord() != null) {
                        for (String bcoord : fusions.Bcoord()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBCOORD,
                                    MOLECULARMATCHMUTATIONSFUSIONBCOORD.BCOORD,
                                    MOLECULARMATCHMUTATIONSFUSIONBCOORD.MOLECULARMATCHMUTATIONSFUSIONID)
                                    .values(bcoord, idFusions)
                                    .execute();
                        }
                    }

                    if (fusions.Bgreg() != null) {
                        for (MolecularMatchBreg breg : fusions.Bgreg()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBREG,
                                    MOLECULARMATCHMUTATIONSFUSIONBREG.NUM,
                                    MOLECULARMATCHMUTATIONSFUSIONBREG.TYPE,
                                    MOLECULARMATCHMUTATIONSFUSIONBREG.MOLECULARMATCHMUTATIONSFUSIONID)
                                    .values(breg.num(), breg.type(), idFusions)
                                    .execute();
                        }
                    }

                    if (fusions.Agreg() != null) {
                        for (MolecularMatchAreg areg : fusions.Agreg()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONAREG,
                                    MOLECULARMATCHMUTATIONSFUSIONAREG.NUM,
                                    MOLECULARMATCHMUTATIONSFUSIONAREG.TYPE,
                                    MOLECULARMATCHMUTATIONSFUSIONAREG.MOLECULARMATCHMUTATIONSFUSIONID)
                                    .values(areg.num(), areg.type(), idFusions)
                                    .execute();
                        }
                    }
                }
            }

            MolecularMatchExonInfo info = mutations.exonsInfo();
            if (info != null) {
                int idExons = context.insertInto(MOLECULARMATCHMUTATIONSEXONSINFO,
                        MOLECULARMATCHMUTATIONSEXONSINFO.TXSTART,
                        MOLECULARMATCHMUTATIONSEXONSINFO.CDSEND,
                        MOLECULARMATCHMUTATIONSEXONSINFO.CHR,
                        MOLECULARMATCHMUTATIONSEXONSINFO.CDSSTART,
                        MOLECULARMATCHMUTATIONSEXONSINFO.TRANSCRIPT,
                        MOLECULARMATCHMUTATIONSEXONSINFO.TXEND,
                        MOLECULARMATCHMUTATIONSEXONSINFO.MOLECULARMATCHMUTATIONSID)
                        .values(info.txStart(), info.cdsEnd(), info.chr(), info.cdsStart(), info.transcript(), info.txEnd(), idMutations)
                        .returning(MOLECULARMATCHMUTATIONSEXONSINFO.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSEXONSINFO.ID);

                int idBoundries = context.insertInto(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES,
                        MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES.MOLECULARMATCHMUTATIONSEXONSINFOID)
                        .values(idExons)
                        .returning(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES.ID);

                int idExonPosities = context.insertInto(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES,
                        MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES.MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESID)
                        .values(idBoundries)
                        .returning(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES.ID);

                if (info.exonBoundaries().exon1() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon1().start(), info.exonBoundaries().exon1().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon2() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon2().start(), info.exonBoundaries().exon2().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon3() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon3().start(), info.exonBoundaries().exon3().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon4() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon4().start(), info.exonBoundaries().exon4().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon5() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon5().start(), info.exonBoundaries().exon5().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon6() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon6().start(), info.exonBoundaries().exon6().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon7() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon7().start(), info.exonBoundaries().exon7().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon8() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon8().start(), info.exonBoundaries().exon8().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon9() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon9().start(), info.exonBoundaries().exon9().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon10() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon10().start(), info.exonBoundaries().exon10().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon11() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon11().start(), info.exonBoundaries().exon11().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon12() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon12().start(), info.exonBoundaries().exon12().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon13() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon13().start(), info.exonBoundaries().exon13().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon14() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon14().start(), info.exonBoundaries().exon14().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon15() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon15().start(), info.exonBoundaries().exon15().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon16() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon16().start(), info.exonBoundaries().exon16().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon17() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon17().start(), info.exonBoundaries().exon17().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon18() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon18().start(), info.exonBoundaries().exon18().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon19() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon19().start(), info.exonBoundaries().exon19().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon20() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon20().start(), info.exonBoundaries().exon20().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon21() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon21().start(), info.exonBoundaries().exon21().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon22() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon22().start(), info.exonBoundaries().exon22().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon23() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon23().start(), info.exonBoundaries().exon23().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon24() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon24().start(), info.exonBoundaries().exon24().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon25() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon25().start(), info.exonBoundaries().exon25().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon26() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon26().start(), info.exonBoundaries().exon26().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon27() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon27().start(), info.exonBoundaries().exon27().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon28() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon28().start(), info.exonBoundaries().exon28().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon29() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon29().start(), info.exonBoundaries().exon29().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon30() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon30().start(), info.exonBoundaries().exon30().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon31() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon31().start(), info.exonBoundaries().exon31().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon32() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon32().start(), info.exonBoundaries().exon32().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon33() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon33().start(), info.exonBoundaries().exon33().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon34() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon34().start(), info.exonBoundaries().exon34().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon35() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon35().start(), info.exonBoundaries().exon35().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon36() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon36().start(), info.exonBoundaries().exon36().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon37() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon37().start(), info.exonBoundaries().exon37().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon38() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon38().start(), info.exonBoundaries().exon38().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon39() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon39().start(), info.exonBoundaries().exon39().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon40() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon40().start(), info.exonBoundaries().exon40().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon41() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon41().start(), info.exonBoundaries().exon41().stop(), idExonPosities)
                            .execute();
                }

            }
        }

    }

    private void writeAssociation(int viccEntryId, @NotNull Association association) {
        int id = context.insertInto(ASSOCIATION,
                ASSOCIATION.VARIANTNAME,
                ASSOCIATION.EVIDENCELEVEL,
                ASSOCIATION.EVIDENCELABEL,
                ASSOCIATION.RESPONSETYPE,
                ASSOCIATION.DRUGLABEL,
                ASSOCIATION.SOURCELINK,
                ASSOCIATION.DESCRIPTION,
                ASSOCIATION.ONCOGENIC,
                ASSOCIATION.VICCENTRYID)
                .values(association.variantName(),
                        association.evidenceLevel(),
                        association.evidenceLabel(),
                        association.responseType(),
                        association.drugLabels(),
                        association.sourceLink(),
                        association.description(),
                        association.oncogenic(),
                        viccEntryId)
                .returning(ASSOCIATION.ID)
                .fetchOne()
                .getValue(ASSOCIATION.ID);
        writeEvidence(id, association.evidence());
        writePublicationsUrls(id, association.publicationUrls());
        writePhenotype(id, association.phenotype());
        writeEnvironmentalContexts(id, association.environmentalContexts());
    }

    private void writeEnvironmentalContexts(int associationId, @Nullable List<EnvironmentalContext> environmentalContexts) {
        if (environmentalContexts != null) {
            for (EnvironmentalContext environmentalContext : environmentalContexts) {
                int id = context.insertInto(ENVIRONMENTALCONTEXT,
                        ENVIRONMENTALCONTEXT.TERM,
                        ENVIRONMENTALCONTEXT.DESCRIPTION,
                        ENVIRONMENTALCONTEXT.SOURCE,
                        ENVIRONMENTALCONTEXT.USANSTEM,
                        ENVIRONMENTALCONTEXT.IDENVIRONMENTALCONTEXT,
                        ENVIRONMENTALCONTEXT.ASSOCIATIONID)
                        .values(environmentalContext.term(),
                                environmentalContext.description(),
                                environmentalContext.source(),
                                environmentalContext.usanStem(),
                                environmentalContext.id(),
                                associationId)
                        .returning(ENVIRONMENTALCONTEXT.ID)
                        .fetchOne()
                        .getValue(ENVIRONMENTALCONTEXT.ID);
                writeApprovedCountries(id, environmentalContext.approvedCountries());
                writeTaxonomy(id, environmentalContext.taxonomy());
            }
        }
    }

    private void writeTaxonomy(int environmentalContextsId, @Nullable Taxonomy taxonomy) {
        if (taxonomy != null) {
            context.insertInto(TAXONOMY,
                    TAXONOMY.KINGDOM,
                    TAXONOMY.DIRECTPARENT,
                    TAXONOMY.CLASS,
                    TAXONOMY.SUBCLASS,
                    TAXONOMY.SUPERCLASS,
                    TAXONOMY.ENVIRONMENTCONTEXTID)
                    .values(taxonomy.kingdom(),
                            taxonomy.directParent(),
                            taxonomy.classs(),
                            taxonomy.subClass(),
                            taxonomy.superClass(),
                            environmentalContextsId)
                    .execute();
        }
    }

    private void writeApprovedCountries(int environmentalContextsId, @NotNull List<String> approvedCountries) {
        for (String approvesCountry : approvedCountries) {
            context.insertInto(APPROVEDCOUNTRY, APPROVEDCOUNTRY.APPROVEDCOUNTRYNAME, APPROVEDCOUNTRY.ENVIRONMENTCONTEXTID)
                    .values(approvesCountry, environmentalContextsId)
                    .execute();
        }
    }

    private void writePhenotype(int associationId, @Nullable Phenotype phenotype) {
        if (phenotype != null) {
            int id = context.insertInto(PHENOTYPE, PHENOTYPE.DESCRIPTION, PHENOTYPE.FAMILY, PHENOTYPE.IDPHENOTYPE, PHENOTYPE.ASSOCIATIONID)
                    .values(phenotype.description(), phenotype.family(), phenotype.id(), associationId)
                    .returning(PHENOTYPE.ID)
                    .fetchOne()
                    .getValue(PHENOTYPE.ID);
            writePhenotypeType(id, phenotype.type());
        }
    }

    private void writePhenotypeType(int phenotypeId, @Nullable PhenotypeType phenotypeType) {
        if (phenotypeType != null) {
            context.insertInto(PHENOTYPETYPE,
                    PHENOTYPETYPE.SOURCE,
                    PHENOTYPETYPE.TERM,
                    PHENOTYPETYPE.IDPHENOTYPETYPE,
                    PHENOTYPETYPE.PHENOTYPEID)
                    .values(phenotypeType.source(), phenotypeType.term(), phenotypeType.id(), phenotypeId)
                    .execute();
        }
    }

    private void writePublicationsUrls(int associationId, @Nullable List<String> publicationsUrls) {
        if (publicationsUrls != null) {
            for (String publicationUrl : publicationsUrls) {
                context.insertInto(PUBLICATIONURL, PUBLICATIONURL.URLOFPUBLICATION, PUBLICATIONURL.ASSOCIATIONID)
                        .values(publicationUrl, associationId)
                        .execute();
            }
        }
    }

    private void writeEvidence(int associationId, @NotNull List<Evidence> evidences) {
        for (Evidence evidence : evidences) {
            int id = context.insertInto(EVIDENCE, EVIDENCE.DESCRIPTION, EVIDENCE.ASSOCIATIONID)
                    .values(evidence.description(), associationId)
                    .returning(EVIDENCE.ID)
                    .fetchOne()
                    .getValue(EVIDENCE.ID);
            writeEvidenceInfo(id, evidence.info());
            writeEvidenceType(id, evidence.evidenceType());
        }
    }

    private void writeEvidenceInfo(int evidenceId, @Nullable EvidenceInfo evidenceInfo) {
        if (evidenceInfo != null) {
            for (String publication : evidenceInfo.publications()) {
                context.insertInto(EVIDENCEINFO, EVIDENCEINFO.PUBLICATION, EVIDENCEINFO.EVIDENCEID)
                        .values(publication, evidenceId)
                        .execute();
            }
        }
    }

    private void writeEvidenceType(int evidenceId, @NotNull EvidenceType evidenceType) {
        context.insertInto(EVIDENCETYPE, EVIDENCETYPE.SOURCENAME, EVIDENCETYPE.IDEVIDENCETYPE, EVIDENCETYPE.EVIDENCEID)
                .values(evidenceType.sourceName(), evidenceType.id(), evidenceId)
                .execute();
    }

    private void writeFeatures(int viccEntryId, @NotNull List<Feature> features) {
        for (Feature feature : features) {
            int id = context.insertInto(FEATURE,
                    FEATURE.NAME,
                    FEATURE.BIOMARKERTYPE,
                    FEATURE.REFERENCENAME,
                    FEATURE.CHROMOSOME,
                    FEATURE.START,
                    FEATURE.END,
                    FEATURE.REF,
                    FEATURE.ALT,
                    FEATURE.PROVENANCERULE,
                    FEATURE.GENESYMBOL,
                    FEATURE.ENTREZID,
                    FEATURE.DESCRIPTION,
                    FEATURE.VICCENTRYID)
                    .values(feature.name(),
                            feature.biomarkerType(),
                            feature.referenceName(),
                            feature.chromosome(),
                            feature.start(),
                            feature.end(),
                            feature.ref(),
                            feature.alt(),
                            feature.provenanceRule(),
                            feature.geneSymbol(),
                            feature.entrezId(),
                            feature.description(),
                            viccEntryId)
                    .returning(FEATURE.ID)
                    .fetchOne()
                    .getValue(FEATURE.ID);
            writeProvenance(id, feature.provenance());
            writeSynonyms(id, feature.synonyms());
            writeLinks(id, feature.links());
            writeSequenceOntology(id, feature.sequenceOntology());
        }
    }

    private void writeProvenance(int featureId, @NotNull List<String> provenances) {
        for (String provenance : provenances) {
            context.insertInto(PROVENANCE, PROVENANCE.PROVENANCENAME, PROVENANCE.FEATUREID).values(provenance, featureId).execute();
        }
    }

    private void writeSynonyms(int featureId, @Nullable List<String> synonyms) {
        if (synonyms != null) {
            for (String synonym : synonyms) {
                context.insertInto(SYNONYM, SYNONYM.SYNONYMNAME, SYNONYM.FEATUREID).values(synonym, featureId).execute();
            }
        }
    }

    private void writeLinks(int featureId, @Nullable List<String> links) {
        if (links != null) {
            for (String link : links) {
                context.insertInto(LINK, LINK.LINKNAME, LINK.FEATUREID).values(link, featureId).execute();
            }
        }
    }

    private void writeSequenceOntology(int featureId, @Nullable SequenceOntology sequenceOntologies) {
        if (sequenceOntologies != null) {
            int id = context.insertInto(SEQUENCEONTOLOGY,
                    SEQUENCEONTOLOGY.SOID,
                    SEQUENCEONTOLOGY.PARENTSOID,
                    SEQUENCEONTOLOGY.NAME,
                    SEQUENCEONTOLOGY.PARENTNAME,
                    SEQUENCEONTOLOGY.FEATUREID)
                    .values(sequenceOntologies.soid(),
                            sequenceOntologies.parentSoid(),
                            sequenceOntologies.name(),
                            sequenceOntologies.parentName(),
                            featureId)
                    .returning(SEQUENCEONTOLOGY.ID)
                    .fetchOne()
                    .getValue(SEQUENCEONTOLOGY.ID);
            writeHierarchy(id, sequenceOntologies.hierarchy());
        }
    }

    private void writeHierarchy(int sequenceOntologyId, @Nullable List<String> hierarchies) {
        if (hierarchies != null) {
            for (String hierarchy : hierarchies) {
                context.insertInto(HIERARCHY, HIERARCHY.HIERARCHYNAME, HIERARCHY.SEQUENCEONTOLOGYID)
                        .values(hierarchy, sequenceOntologyId)
                        .execute();
            }
        }
    }

    private void writeTags(int viccEntryId, @NotNull List<String> tags) {
        for (String tag : tags) {
            context.insertInto(TAG, TAG.TAGNAME, TAG.VICCENTRYID).values(tag, viccEntryId).execute();
        }
    }

    private void writeDevTags(int viccEntryId, @NotNull List<String> devTags) {
        for (String devTag : devTags) {
            context.insertInto(DEVTAG, DEVTAG.DEVTAGNAME, DEVTAG.VICCENTRYID).values(devTag, viccEntryId).execute();
        }
    }

    private void writeGeneIdentifiers(int viccEntryId, @NotNull List<GeneIdentifier> geneIdentifiers) {
        for (GeneIdentifier geneIdentifier : geneIdentifiers) {
            context.insertInto(GENEIDENTIFIER,
                    GENEIDENTIFIER.SYMBOL,
                    GENEIDENTIFIER.ENTREZID,
                    GENEIDENTIFIER.ENSEMBLGENEID,
                    GENEIDENTIFIER.VICCENTRYID)
                    .values(geneIdentifier.symbol(), geneIdentifier.entrezId(), geneIdentifier.ensemblGeneId(), viccEntryId)
                    .execute();
        }
    }

    private void writeGenes(int viccEntryId, @NotNull List<String> genes) {
        for (String gene : genes) {
            context.insertInto(GENE, GENE.GENENAME, GENE.VICCENTRYID).values(gene, viccEntryId).execute();
        }
    }

    private void writeFeatureNames(int viccEntryId, @Nullable List<String> featureNames) {
        if (featureNames != null) {
            for (String featureName : featureNames) {
                context.insertInto(FEATURENAME, FEATURENAME.NAMEOFFEATURE, FEATURENAME.VICCENTRYID)
                        .values(featureName, viccEntryId)
                        .execute();
            }
        }
    }

    public void deleteAll() {
        LOGGER.info("Deleting all from vicc db");

        context.deleteFrom(TAG).execute();
        context.deleteFrom(DEVTAG).execute();
        context.deleteFrom(GENEIDENTIFIER).execute();
        context.deleteFrom(GENE).execute();
        context.deleteFrom(FEATURENAME).execute();
        context.deleteFrom(FEATURE).execute();
        context.deleteFrom(PROVENANCE).execute();
        context.deleteFrom(SYNONYM).execute();
        context.deleteFrom(LINK).execute();
        context.deleteFrom(SEQUENCEONTOLOGY).execute();
        context.deleteFrom(HIERARCHY).execute();
        context.deleteFrom(ASSOCIATION).execute();
        context.deleteFrom(EVIDENCE).execute();
        context.deleteFrom(EVIDENCEINFO).execute();
        context.deleteFrom(EVIDENCETYPE).execute();
        context.deleteFrom(PUBLICATIONURL).execute();
        context.deleteFrom(PHENOTYPE).execute();
        context.deleteFrom(PHENOTYPETYPE).execute();
        context.deleteFrom(ENVIRONMENTALCONTEXT).execute();
        context.deleteFrom(APPROVEDCOUNTRY).execute();
        context.deleteFrom(TAXONOMY).execute();
        context.deleteFrom(VICCENTRY).execute();
        context.deleteFrom(SAGE).execute();
        context.deleteFrom(BRCAPART1).execute();
        context.deleteFrom(BRCAPART2).execute();
        context.deleteFrom(CGI).execute();
        context.deleteFrom(CGICDNA).execute();
        context.deleteFrom(CGIINDIVIDUALMUTATION).execute();
        context.deleteFrom(CGIGDNA).execute();
        context.deleteFrom(CGITRANSCRIPT).execute();
        context.deleteFrom(CGISTRAND).execute();
        context.deleteFrom(CGIINFO).execute();
        context.deleteFrom(CGIREGION).execute();
        context.deleteFrom(JAX).execute();
        context.deleteFrom(JAXMOLECULARPROFILE).execute();
        context.deleteFrom(JAXTHERAPY).execute();
        context.deleteFrom(JAXINDICATIONS).execute();
        context.deleteFrom(JAXREFERENCES).execute();
        context.deleteFrom(JAXTRIALS).execute();
        context.deleteFrom(JAXTRIALSINDICATIONS).execute();
        context.deleteFrom(JAXTRIALSVARIANTREQUIREMENTDETAILS).execute();
        context.deleteFrom(JAXTRIALSMOLECULARPROFILE).execute();
        context.deleteFrom(JAXTRIALSTHERAPIES).execute();
        context.deleteFrom(PMKB).execute();
        context.deleteFrom(PMKBTISSUE).execute();
        context.deleteFrom(PMKBTUMOR).execute();
        context.deleteFrom(PMKBVARIANT).execute();
        context.deleteFrom(PMKBGENE).execute();
        context.deleteFrom(ONCOKB).execute();
        context.deleteFrom(ONCOKBBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBVARIANTBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBCONSEQUENCESBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBGENEBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBGENEALIASESBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBCLINICAL).execute();
        context.deleteFrom(ONCOKBDRUGABSTRACTSCLINICAL).execute();
        context.deleteFrom(ONCOKBVARIANTCLINICAL).execute();
        context.deleteFrom(ONCOKBCONSEQUENCESCLINICAL).execute();
        context.deleteFrom(ONCOKBGENECLINICAL).execute();
        context.deleteFrom(ONCOKBGENEALIASESCLINICAL).execute();
        context.deleteFrom(MOLECULARMATCHTRIALS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSALTERATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSINTERVATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOTHERNAME).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOTHERGROUPLABEL).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOVERALLCONTACT).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSTAGS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSLOCATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSCONTACT).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSGEO).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSLOCATION).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSCOORDINATES).execute();
        context.deleteFrom(CIVIC).execute();
        context.deleteFrom(CIVICASSERTIONS).execute();
        context.deleteFrom(CIVICHGVSEXPRESSIONS).execute();
        context.deleteFrom(CIVICCLINVARENTRIES).execute();
        context.deleteFrom(CIVICVARIANTALIASES).execute();
        context.deleteFrom(CIVICVARIANTTYPES).execute();
        context.deleteFrom(CIVICDESCRIPTION).execute();
        context.deleteFrom(CIVICCOORDINATES).execute();
        context.deleteFrom(CIVICVARIANTSGROUPS).execute();
        context.deleteFrom(CIVICVARIANTSGROUPSCOORDINATES).execute();
        context.deleteFrom(CIVICVARIANTSGROUPSTYPES).execute();
        context.deleteFrom(CIVICEVIDENCEITEMS).execute();
        context.deleteFrom(CIVICDRUGS).execute();
        context.deleteFrom(CIVICDISEASE).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSSOURCE).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSPUBLICATION).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSCLINICALTRIAL).execute();
        context.deleteFrom(CIVICSOURCE).execute();
        context.deleteFrom(CIVICERROR).execute();
        context.deleteFrom(CIVICPUBLICATION).execute();
        context.deleteFrom(CIVICCLINICALTRIAL).execute();
        context.deleteFrom(CIVICLIFECYCLEACTIONS).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDON).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONUSER).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONAVATARS).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONORGANIZATION).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTMODIFIED).execute();
        context.deleteFrom(CIVICLASTMODIFIEDUSER).execute();
        context.deleteFrom(CIVICLASTMODIFIEDAVATARS).execute();
        context.deleteFrom(CIVICLASTMODIFIEDORGANIZATION).execute();
        context.deleteFrom(CIVICLASTMODIFIEDPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTREVIEWED).execute();
        context.deleteFrom(CIVICLASTREVIEWEDUSER).execute();
        context.deleteFrom(CIVICLASTREVIEWEDAVATARS).execute();
        context.deleteFrom(CIVICLASTREVIEWEDORGANIZATION).execute();
        context.deleteFrom(CIVICLASTREVIEWEDPROFILEIMAGE).execute();
        context.deleteFrom(MOLECULARMATCH).execute();
        context.deleteFrom(MOLECULARMATCHINSTUTITION).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEGENE0).execute();
        context.deleteFrom(MOLECULARMATCHEXTERNALID).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDESTAGE0).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEDRUG1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDECONDITION1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEMUTATION1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDECONDITION0).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEMUTATION0).execute();
        context.deleteFrom(MOLECULARMATCHCRITERIAMET).execute();
        context.deleteFrom(MOLECULARMATCHSOURCE).execute();
        context.deleteFrom(MOLECULARMATCHTIEREXPLANATION).execute();
        context.deleteFrom(MOLECULARMATCHTHERAPEUTICCONTEXT).execute();
        context.deleteFrom(MOLECULARMATCHTAGS).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFO).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOCONSEQUENCES).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOFUSIONS).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOLOCATIONS).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATION).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEND).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONSTART).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONCHR).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPATHOLOGY).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONREF).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEXON).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONALT).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONSOURCES).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONCOSMICID).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONDBSNP).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPOPFREQMAX).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEXONICFUNC).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPARENTS).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS).execute();
        context.deleteFrom(MOLECULARMATCHCRITERIAUNMET).execute();
        context.deleteFrom(MOLECULARMATCHPREFELANCE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSMUTATIONTYPE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSSOURCE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSSYNONYMS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSPATHOLOGY).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSCDNA).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSPARENTS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATA).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATAFULLAA).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATAGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSMAP).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSGRCH37LOCATION).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSION).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBCHR).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONAGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBTX).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONACHR).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONINS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONACOORD).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBORI).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONABAND).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBBAND).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONAORI).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONATX).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBCOORD).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBREG).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONAREG).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONSINFO).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41).execute();
        //TODO: adding object asts of molecular match
    }
}
