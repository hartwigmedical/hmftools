package com.hartwig.hmftools.patientdb.log4j2plugins;

import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;

import org.apache.logging.log4j.core.LogEvent;
import org.apache.logging.log4j.core.config.plugins.Plugin;
import org.apache.logging.log4j.core.pattern.ConverterKeys;
import org.apache.logging.log4j.core.pattern.LogEventPatternConverter;

@Plugin(name = "PatientIdConverter",
        category = "Converter")
@ConverterKeys({ "PatientId" })
public class PatientIdPatternConverter extends LogEventPatternConverter {
    protected PatientIdPatternConverter(String name, String style) {
        super(name, style);
    }

    public static PatientIdPatternConverter newInstance(String[] options) {
        return new PatientIdPatternConverter("PatientId", "PatientId");
    }

    @Override
    public void format(LogEvent event, StringBuilder toAppendTo) {
        if (event.getMessage() instanceof ValidationFinding) {
            final ValidationFinding dbLogEvent = (ValidationFinding) event.getMessage();
            toAppendTo.append(dbLogEvent.patientId());
        }
    }
}