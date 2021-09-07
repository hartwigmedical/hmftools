package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfDataField;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfResolveException;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatusKey;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatusModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.ImmutableFormStatusKey;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class XMLPatientReader extends EcrfReader {

    private static final Logger LOGGER = LogManager.getLogger(XMLPatientReader.class);

    private static final String PATIENT_TAG = "SubjectData";
    private static final String PATIENT_ID_ATTRIBUTE = "SubjectKey";

    private static final String STUDY_EVENT_TAG = "StudyEventData";
    private static final String STUDY_EVENT_OID_ATTRIBUTE = "StudyEventOID";
    private static final String STUDY_EVENT_REPEAT_KEY_ATTRIBUTE = "StudyEventRepeatKey";
    private static final String FORM_TAG = "FormData";
    private static final String FORM_OID_ATTRIBUTE = "FormOID";
    private static final String FORM_REPEAT_KEY_ATTRIBUTE = "FormRepeatKey";
    private static final String ITEM_GROUP_TAG = "ItemGroupData";
    private static final String ITEM_GROUP_OID_ATTRIBUTE = "ItemGroupOID";
    private static final String ITEM_GROUP_REPEAT_KEY = "ItemGroupRepeatKey";

    private static final String ITEM_TAG = "ItemData";
    private static final String ITEM_OID_ATTRIBUTE = "ItemOID";
    private static final String ITEM_VALUE_ATTRIBUTE = "Value";

    // TODO See DEV-491, DEV-1304, DEV-1724
    private static final List DROPDOWN_ITEMS_TO_IGNORE_FOR_WARNING = Lists.newArrayList("FLD.DDELINV",
            "FLD.DDELINST",
            "FLD.CS.CSRTYPE",
            "FLD.CS.CSRSTYPE",
            "FLD.CARCINOMA.ENTRYSTAGE",
            "FLD.TRTAFTER.PLANNEDTRT",
            "FLD.TRTAFTER.ENDTRTREAS",
            "FLD.TRTAFTER.CPCTTRIALNR",
            "FLD.INFORMEDCONSENT.ADDNA",
            "FLD.INFORMEDCONSENT.ADPID",
            "FLD.INFORMEDCONSENT.ADGEN");

    private XMLPatientReader() {
    }

    @NotNull
    public static List<EcrfPatient> readPatients(@NotNull XMLStreamReader reader, @NotNull final XMLEcrfDatamodel datamodel,
            @NotNull final FormStatusModel formStatusModel) throws XMLStreamException {
        List<EcrfPatient> patients = Lists.newArrayList();

        while (reader.hasNext() && !isClinicalDataEnd(reader)) {
            if (isPatientStart(reader)) {
                patients.add(readPatient(reader, datamodel, formStatusModel));
            }
            reader.next();
        }

        return patients;
    }

    @NotNull
    private static EcrfPatient readPatient(@NotNull final XMLStreamReader reader, @NotNull final XMLEcrfDatamodel datamodel,
            @NotNull final FormStatusModel formStatusModel) throws XMLStreamException {
        String patientId = reformatPatientId(reader.getAttributeValue("", PATIENT_ID_ATTRIBUTE));
        Map<String, List<EcrfStudyEvent>> studyEventsPerOID = Maps.newHashMap();
        List<EcrfDataField> fields = Lists.newArrayList();

        String currentStudyEventOID = Strings.EMPTY;
        String currentStudyEventIdx = Strings.EMPTY;
        String currentFormOID = Strings.EMPTY;
        String currentFormIdx = Strings.EMPTY;
        String currentItemGroupOID = Strings.EMPTY;
        String currentItemGroupIdx = Strings.EMPTY;
        EcrfStudyEvent currentStudy = new EcrfStudyEvent();
        EcrfForm currentForm = new EcrfForm(FormStatus.undefined());
        EcrfItemGroup currentItemGroup = new EcrfItemGroup();
        while (reader.hasNext() && !isPatientEnd(reader)) {
            if (isStudyEventStart(reader)) {
                currentStudyEventOID = reader.getAttributeValue("", STUDY_EVENT_OID_ATTRIBUTE);
                currentStudyEventIdx = reader.getAttributeValue("", STUDY_EVENT_REPEAT_KEY_ATTRIBUTE);
                currentStudy = new EcrfStudyEvent();
                if (!studyEventsPerOID.containsKey(currentStudyEventOID)) {
                    studyEventsPerOID.put(currentStudyEventOID, Lists.newArrayList());
                }
                studyEventsPerOID.get(currentStudyEventOID).add(currentStudy);
            } else if (isFormStart(reader)) {
                currentFormOID = reader.getAttributeValue("", FORM_OID_ATTRIBUTE);
                currentFormIdx = reader.getAttributeValue("", FORM_REPEAT_KEY_ATTRIBUTE);
                final String formName = datamodel.forms().get(currentFormOID).name();
                final String studyEventName = datamodel.studyEvents().get(currentStudyEventOID).name();
                final FormStatusKey formKey =
                        new ImmutableFormStatusKey(patientId, formName, currentFormIdx, studyEventName, currentStudyEventIdx);
                final FormStatus formStatus = formStatusModel.formStatuses().get(formKey);
                currentForm = formStatus != null ? new EcrfForm(formStatus) : new EcrfForm(FormStatus.undefined());
                currentStudy.addForm(currentFormOID, currentForm);
            } else if (isItemGroupStart(reader)) {
                currentItemGroupOID = reader.getAttributeValue("", ITEM_GROUP_OID_ATTRIBUTE);
                currentItemGroupIdx = reader.getAttributeValue("", ITEM_GROUP_REPEAT_KEY);
                currentItemGroup = new EcrfItemGroup();
                currentForm.addItemGroup(currentItemGroupOID, currentItemGroup);
            } else if (isFieldStart(reader)) {
                String OID = reader.getAttributeValue("", ITEM_OID_ATTRIBUTE);
                String value = Strings.EMPTY;
                try {
                    value = datamodel.resolveValue(OID, reader.getAttributeValue("", ITEM_VALUE_ATTRIBUTE));
                } catch (EcrfResolveException exception) {
                    if (!DROPDOWN_ITEMS_TO_IGNORE_FOR_WARNING.contains(OID)) {
                        LOGGER.warn("Resolve issue for {}: {}", patientId, exception.getMessage());
                    }
                }
                currentItemGroup.addItem(OID, value);
                fields.add(EcrfDataField.of(patientId,
                        currentStudyEventOID,
                        currentStudyEventIdx,
                        currentFormOID,
                        currentFormIdx,
                        currentItemGroupOID,
                        currentItemGroupIdx,
                        OID,
                        value,
                        currentForm.status().stateString(),
                        Boolean.toString(currentForm.locked())));
            }
            reader.next();
        }
        return new EcrfPatient(patientId, studyEventsPerOID, fields);
    }

    @NotNull
    private static String reformatPatientId(@NotNull String ecrfPatientId) {
        return ecrfPatientId.replaceAll("-", "");
    }

    private static boolean isPatientStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, PATIENT_TAG);
    }

    private static boolean isPatientEnd(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, PATIENT_TAG);
    }

    private static boolean isFieldStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_TAG);
    }

    private static boolean isStudyEventStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, STUDY_EVENT_TAG);
    }

    private static boolean isFormStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, FORM_TAG);
    }

    private static boolean isItemGroupStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_GROUP_TAG);
    }
}
