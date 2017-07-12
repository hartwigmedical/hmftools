package com.hartwig.hmftools.patientdb.log4j2plugins;

import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;

import org.apache.logging.log4j.core.LogEvent;
import org.apache.logging.log4j.core.config.plugins.Plugin;
import org.apache.logging.log4j.core.pattern.ConverterKeys;
import org.apache.logging.log4j.core.pattern.LogEventPatternConverter;

@Plugin(name = "FormStatusConverter",
        category = "Converter")
@ConverterKeys({ "FormStatus" })
public class FormStatusPatternConverter extends LogEventPatternConverter {
    protected FormStatusPatternConverter(String name, String style) {
        super(name, style);
    }

    public static FormStatusPatternConverter newInstance(String[] options) {
        return new FormStatusPatternConverter("FormStatus", "FormStatus");
    }

    @Override
    public void format(LogEvent event, StringBuilder toAppendTo) {
        if (event.getMessage() instanceof ValidationFinding) {
            final ValidationFinding dbLogEvent = (ValidationFinding) event.getMessage();
            toAppendTo.append(dbLogEvent.formStatus());
        }
    }
}
