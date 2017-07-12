package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.ECRF;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.EcrfItem;
import com.hartwig.hmftools.patientdb.data.ImmutableEcrfItem;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class EcrfDAO {
    private static final Logger LOGGER = LogManager.getLogger(EcrfDAO.class);

    @NotNull
    private final DSLContext context;

    EcrfDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void clear() {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(ECRF).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    void writePatient(@NotNull final EcrfPatient patient, final boolean sequenced) {
        LOGGER.info("writing patient: " + patient.patientId());
        final List<EcrfItem> ecrfItems = patientToEcrfItems(patient, sequenced);
        context.batch(ecrfItems.stream()
                .map(ecrfItem -> context.insertInto(ECRF, ECRF.PATIENTID, ECRF.STUDYEVENT, ECRF.STUDYEVENTIDX, ECRF.FORM, ECRF.FORMIDX,
                        ECRF.ITEMGROUP, ECRF.ITEMGROUPIDX, ECRF.ITEM, ECRF.ITEMVALUE, ECRF.FORMSTATUS, ECRF.LOCKED, ECRF.SEQUENCED)
                        .values(ecrfItem.patientId(), ecrfItem.studyOID(), ecrfItem.studyIdx(), ecrfItem.formOID(), ecrfItem.formIdx(),
                                ecrfItem.itemGroupOID(), ecrfItem.itemGroupIdx(), ecrfItem.itemOID(), ecrfItem.itemValue(),
                                ecrfItem.formStatus(), ecrfItem.locked(), ecrfItem.sequencedString()))
                .collect(Collectors.toList())).execute();
        LOGGER.info("finished patient: " + patient.patientId());
    }

    @NotNull
    private List<EcrfItem> patientToEcrfItems(@NotNull final EcrfPatient patient, final boolean sequenced) {
        final String patientId = patient.patientId();
        final List<EcrfItem> ecrfItems = Lists.newArrayList();
        for (final String studyOID : patient.studyEventsPerOID().keySet()) {
            final List<EcrfStudyEvent> studyEvents = patient.studyEventsPerOID(studyOID);
            for (int studyIdx = 0; studyIdx < studyEvents.size(); studyIdx++) {
                for (final String formOID : studyEvents.get(studyIdx).formsPerOID().keySet()) {
                    final List<EcrfForm> forms = studyEvents.get(studyIdx).formsPerOID().get(formOID);
                    for (int formIdx = 0; formIdx < forms.size(); formIdx++) {
                        final String formStatus = forms.get(formIdx).status();
                        final String formLocked = forms.get(formIdx).locked();
                        for (final String itemGroupOID : forms.get(formIdx).itemGroupsPerOID().keySet()) {
                            final List<EcrfItemGroup> itemGroups = forms.get(formIdx).itemGroupsPerOID().get(itemGroupOID);
                            for (int itemGroupIdx = 0; itemGroupIdx < itemGroups.size(); itemGroupIdx++) {
                                for (final String itemOID : itemGroups.get(itemGroupIdx).itemsPerOID().keySet()) {
                                    final List<String> items = itemGroups.get(itemGroupIdx).itemsPerOID().get(itemOID);
                                    if (items.size() > 0) {
                                        final String itemValue = items.get(0);
                                        for (final String item : items) {
                                            if (!itemValue.equals(item)) {
                                                LOGGER.warn(
                                                        "Not all items in item group: " + itemGroupOID + " equal: " + item + "; expected: "
                                                                + itemValue);
                                                break;
                                            }
                                        }
                                        ecrfItems.add(new ImmutableEcrfItem(patientId, studyOID, studyIdx, formOID, formIdx, itemGroupOID,
                                                itemGroupIdx, itemOID, itemValue, sequenced, formStatus, formLocked));

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return ecrfItems;
    }
}
