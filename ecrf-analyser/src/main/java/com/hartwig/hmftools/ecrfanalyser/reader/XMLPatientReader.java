package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfFieldFunctions;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfPatient;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfResolveException;

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
    private static final String FIELD_TAG = "ItemData";
    private static final String FIELD_OID_ATTRIBUTE = "ItemOID";
    private static final String FIELD_VALUE_ATTRIBUTE = "Value";

    private XMLPatientReader() {
    }

    @NotNull
    public static List<EcrfPatient> readPatients(@NotNull XMLStreamReader reader, @NotNull Iterable<EcrfField> fields)
            throws XMLStreamException {
        List<EcrfPatient> patients = Lists.newArrayList();
        Map<String, EcrfField> OIDToEcrfFields = mapOIDToEcrfFields(fields);

        while (reader.hasNext() && !isClinicalDataEnd(reader)) {
            if (isPatientStart(reader)) {
                patients.add(readPatient(reader, OIDToEcrfFields));
            }
            reader.next();
        }

        return patients;
    }

    @NotNull
    private static Map<String, EcrfField> mapOIDToEcrfFields(@NotNull Iterable<EcrfField> fields) {
        Map<String, EcrfField> mapping = Maps.newHashMap();
        for (EcrfField field : fields) {
            mapping.put(field.itemOID(), field);
        }
        return mapping;
    }

    @NotNull
    private static EcrfPatient readPatient(@NotNull XMLStreamReader reader,
            @NotNull Map<String, EcrfField> OIDtoEcrfFieldMap) throws XMLStreamException {
        String patientId = toCPCTPatientId(reader.getAttributeValue("", PATIENT_ID_ATTRIBUTE));
        Map<EcrfField, String> fieldValues = Maps.newHashMap();

        String currentStudyEvent;
        String currentForm;
        String currentItemGroup;
        while (reader.hasNext() && !isPatientEnd(reader)) {
            if (isStudyEventStart(reader)) {
                currentStudyEvent = reader.getAttributeValue("", STUDY_EVENT_OID_ATTRIBUTE);
            } else if (isFormStart(reader)) {
                currentForm = reader.getAttributeValue("", FORM_OID_ATTRIBUTE);
            } else if (isItemGroupStart(reader)) {
                currentItemGroup = reader.getAttributeValue("", ITEM_GROUP_OID_ATTRIBUTE);
            } else if (isFieldStart(reader)) {
                String OID = reader.getAttributeValue("", FIELD_OID_ATTRIBUTE);
                EcrfField field = OIDtoEcrfFieldMap.get(OID);
                if (field != null) {
                    String value = Strings.EMPTY;
                    try {
                        value = EcrfFieldFunctions.resolveValue(field,
                                reader.getAttributeValue("", FIELD_VALUE_ATTRIBUTE));
                    } catch (EcrfResolveException exception) {
                        LOGGER.warn("Resolve issue for " + patientId + ": " + exception.getMessage());
                    }
                    fieldValues.put(field, value);
                }
            }
            reader.next();
        }
        return new EcrfPatient(patientId, fieldValues);
    }

    @NotNull
    private static String toCPCTPatientId(@NotNull String ecrfPatientId) {
        return ecrfPatientId.replaceAll("-", "");
    }

    private static boolean isPatientStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, PATIENT_TAG);
    }

    private static boolean isPatientEnd(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, PATIENT_TAG);
    }

    private static boolean isFieldStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, FIELD_TAG);
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
