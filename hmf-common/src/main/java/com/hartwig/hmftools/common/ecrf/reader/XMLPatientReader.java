package com.hartwig.hmftools.common.ecrf.reader;

import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfFieldFunctions;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfResolveException;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;

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
    private static final String FORM_TAG = "FormData";
    private static final String FORM_OID_ATTRIBUTE = "FormOID";
    private static final String ITEM_GROUP_TAG = "ItemGroupData";
    private static final String ITEM_GROUP_OID_ATTRIBUTE = "ItemGroupOID";
    private static final String ITEM_TAG = "ItemData";
    private static final String ITEM_OID_ATTRIBUTE = "ItemOID";
    private static final String ITEM_VALUE_ATTRIBUTE = "Value";

    private XMLPatientReader() {
    }

    @NotNull
    public static List<EcrfPatient> readPatients(@NotNull XMLStreamReader reader, @NotNull Iterable<EcrfField> fields)
            throws XMLStreamException {
        final List<EcrfPatient> patients = Lists.newArrayList();
        final Map<String, EcrfField> nameToEcrfFields = mapNameToEcrfFields(fields);

        while (reader.hasNext() && !isClinicalDataEnd(reader)) {
            if (isPatientStart(reader)) {
                patients.add(readPatient(reader, nameToEcrfFields));
            }
            reader.next();
        }

        return patients;
    }

    @NotNull
    private static Map<String, EcrfField> mapNameToEcrfFields(@NotNull Iterable<EcrfField> fields) {
        final Map<String, EcrfField> mapping = Maps.newHashMap();
        for (final EcrfField field : fields) {
            mapping.put(field.name(), field);
        }
        return mapping;
    }

    @NotNull
    private static EcrfPatient readPatient(@NotNull XMLStreamReader reader,
            @NotNull Map<String, EcrfField> nameToEcrfFieldMap) throws XMLStreamException {
        final String patientId = toCPCTPatientId(reader.getAttributeValue("", PATIENT_ID_ATTRIBUTE));
        final Map<EcrfField, List<String>> fieldValues = initializeFieldValues(nameToEcrfFieldMap.values());
        final Map<String, List<EcrfStudyEvent>> studyEventsPerOID = Maps.newHashMap();

        String currentStudyEventOID = Strings.EMPTY;
        String currentFormOID = Strings.EMPTY;
        String currentItemGroupOID = Strings.EMPTY;
        EcrfStudyEvent currentStudy = new EcrfStudyEvent();
        EcrfForm currentForm = new EcrfForm();
        EcrfItemGroup currentItemGroup = new EcrfItemGroup();
        while (reader.hasNext() && !isPatientEnd(reader)) {
            if (isStudyEventStart(reader)) {
                currentStudyEventOID = reader.getAttributeValue("", STUDY_EVENT_OID_ATTRIBUTE);
                currentStudy = new EcrfStudyEvent();
                if (!studyEventsPerOID.containsKey(currentStudyEventOID)) {
                    studyEventsPerOID.put(currentStudyEventOID, Lists.newArrayList());
                }
                studyEventsPerOID.get(currentStudyEventOID).add(currentStudy);
            } else if (isFormStart(reader)) {
                currentFormOID = reader.getAttributeValue("", FORM_OID_ATTRIBUTE);
                currentForm = new EcrfForm();
                currentStudy.addForm(currentFormOID, currentForm);
            } else if (isItemGroupStart(reader)) {
                currentItemGroupOID = reader.getAttributeValue("", ITEM_GROUP_OID_ATTRIBUTE);
                currentItemGroup = new EcrfItemGroup();
                currentForm.addItemGroup(currentItemGroupOID, currentItemGroup);
            } else if (isFieldStart(reader)) {
                String OID = reader.getAttributeValue("", ITEM_OID_ATTRIBUTE);
                String name = EcrfFieldFunctions.name(currentStudyEventOID, currentFormOID, currentItemGroupOID, OID);
                if (EcrfFieldFunctions.isRelevant(name)) {
                    final EcrfField field = nameToEcrfFieldMap.get(name);
                    if (field != null) {
                        String value = Strings.EMPTY;
                        try {
                            value = EcrfFieldFunctions.resolveValue(field,
                                    reader.getAttributeValue("", ITEM_VALUE_ATTRIBUTE));
                        } catch (EcrfResolveException exception) {
                            LOGGER.warn("Resolve issue for " + patientId + ": " + exception.getMessage());
                        }
                        currentItemGroup.addItem(OID, value);
                        fieldValues.get(field).add(value);
                    } else {
                        LOGGER.warn("Could not resolve field with name " + name);
                    }
                }
            }
            reader.next();
        }
        return new EcrfPatient(patientId, fieldValues, studyEventsPerOID);
    }

    @NotNull
    private static Map<EcrfField, List<String>> initializeFieldValues(@NotNull Iterable<EcrfField> fields) {
        final Map<EcrfField, List<String>> fieldValues = Maps.newHashMap();
        for (final EcrfField field : fields) {
            fieldValues.put(field, Lists.newArrayList());
        }
        return fieldValues;
    }

    @NotNull
    private static String toCPCTPatientId(@NotNull final String ecrfPatientId) {
        return ecrfPatientId.replaceAll("-", "");
    }

    private static boolean isPatientStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, PATIENT_TAG);
    }

    private static boolean isPatientEnd(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, PATIENT_TAG);
    }

    private static boolean isFieldStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_TAG);
    }

    private static boolean isStudyEventStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, STUDY_EVENT_TAG);
    }

    private static boolean isFormStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, FORM_TAG);
    }

    private static boolean isItemGroupStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_GROUP_TAG);
    }
}
