package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class EcrfPatientReader extends EcrfReader {

    private static final Logger LOGGER = LogManager.getLogger(EcrfPatientReader.class);

    private static final String PATIENT_TAG = "SubjectData";
    private static final String PATIENT_ID_ATTRIBUTE = "SubjectKey";

    private static final String FIELD_TAG = "ItemData";
    private static final String FIELD_OID_ATTRIBUTE = "ItemOID";
    private static final String FIELD_VALUE_ATTRIBUTE = "Value";

    private EcrfPatientReader() {
    }

    @NotNull
    public static List<EcrfPatient> readPatients(@NotNull XMLStreamReader reader, @NotNull List<EcrfField> fields)
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
    private static Map<String, EcrfField> mapOIDToEcrfFields(@NotNull List<EcrfField> fields) {
        Map<String, EcrfField> mapping = Maps.newHashMap();
        for (EcrfField field : fields) {
            mapping.put(OIDFunctions.toOID(field.category(), field.fieldName()), field);
        }
        return mapping;
    }

    @NotNull
    private static EcrfPatient readPatient(@NotNull XMLStreamReader reader,
            @NotNull Map<String, EcrfField> OIDtoEcrfFieldMap) throws XMLStreamException {
        String patientId = toCPCTPatientId(reader.getAttributeValue("", PATIENT_ID_ATTRIBUTE));
        Map<EcrfField, String> fieldValues = Maps.newHashMap();

        while (reader.hasNext() && !isPatientEnd(reader)) {
            if (isFieldStart(reader)) {
                String OID = reader.getAttributeValue("", FIELD_OID_ATTRIBUTE);
                EcrfField field = OIDtoEcrfFieldMap.get(OID);
                if (field != null) {
                    String value = reader.getAttributeValue("", FIELD_VALUE_ATTRIBUTE);
                    if (!field.isFreeText()) {
                        value = field.values().get(Integer.valueOf(value));
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
}
