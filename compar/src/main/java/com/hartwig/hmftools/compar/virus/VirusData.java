package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.compar.common.CategoryType.VIRUS;

import java.util.List;

import com.hartwig.hmftools.common.viridian.VirusReport;
import com.hartwig.hmftools.common.virusinterpreter.AnnotatedVirus;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.virus.VirusComparer.Fields;

public class VirusData extends ComparableItem
{
    public final String Key;
    public final boolean Reported;

    private VirusData(final String key, final boolean reported)
    {
        Key = key;
        Reported = reported;
    }

    static VirusData fromViridian(final VirusReport report, final List<FieldInfo> fields)
    {
        VirusData data = new VirusData(report.virusGroup(), report.isReported());
        data.addStringValue(Fields.Name.toString(), report.virusName(), fields);
        data.addBoolValue(Fields.Reported.toString(), report.isReported(), fields);
        data.addStringValue(Fields.Interpretation.toString(), String.valueOf(report.reportingType()), fields);
        data.addIntValue(Fields.Integrations.toString(), report.integrations(), fields);
        data.addDoubleValue(Fields.MeanCoverage.toString(), report.meanDepth(), fields);
        data.addDoubleValue(Fields.PercentageCovered.toString(), report.coveragePercentage(), fields);

        if(report.copiesPerTumorCell() != null)
            data.addDoubleValue(Fields.CopiesPerTumorCell.toString(), report.copiesPerTumorCell(), fields);

        data.addStringValue(Fields.ClonalityStatus.toString(), report.clonalityStatus().name(), fields);
        return data;
    }

    static VirusData fromVirusInterpreter(final AnnotatedVirus virus, final List<FieldInfo> fields)
    {
        VirusData data = new VirusData(virusGroup(virus), virus.reported());
        data.addStringValue(Fields.Name.toString(), virus.name(), fields);
        data.addBoolValue(Fields.Reported.toString(), virus.reported(), fields);
        data.addStringValue(Fields.Interpretation.toString(), String.valueOf(virus.interpretation()), fields);
        data.addIntValue(Fields.Integrations.toString(), virus.integrations(), fields);
        data.addDoubleValue(Fields.MeanCoverage.toString(), virus.meanCoverage(), fields);
        data.addDoubleValue(Fields.PercentageCovered.toString(), virus.percentageCovered() / 100, fields);
        return data;
    }

    private static String virusGroup(final AnnotatedVirus virus)
    {
        if(virus.interpretation() != null)
        {
            return switch(virus.interpretation())
            {
                case HPV -> "Human papillomavirus";
                case EBV -> "Human gammaherpesvirus 4";
                case HHV8 -> "Human gammaherpesvirus 8";
                case HBV -> "Hepatitis B virus";
                case MCV -> "Merkel cell polyomavirus";
            };
        }
        else
        {
            String lowerCaseName = virus.name().toLowerCase();
            if(lowerCaseName.startsWith("hepatitis c virus"))
            {
                return "Hepatitis C virus";
            }
            else
            {
                return switch(lowerCaseName)
                {
                    case "hepacivirus hominis" -> "Hepatitis C virus";
                    case "betapolyomavirus hominis" -> "BK polyomavirus";
                    case "human t-cell leukemia virus type i" -> "Human T-lymphotropic virus 1";
                    default -> virus.name();
                };
            }
        }
    }

    @Override
    public CategoryType category()
    {
        return VIRUS;
    }

    @Override
    public String key()
    {
        return Key;
    }

    @Override
    public boolean reportable()
    {
        return Reported;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final VirusData otherData = (VirusData) other;
        return Key.equals(otherData.Key);
    }
}
