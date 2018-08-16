hmf_id_hashes has been generated using v1.2 of hmf-id-generator.

The set of 3611 patient ids was queried on August 16th 2018, 12:10:00 with the following:

select distinct patientIdentifier from
(select patientIdentifier from patient
UNION
SELECT concat('CPCT02', CPCTCNT.itemValue, LPAD(RIGHT(CPCTPN.itemValue,4), 4, '0')) as patientIdentifier
 FROM drupEcrf CTCT2YN,  drupEcrf CPCTCNT, drupEcrf CPCTPN
WHERE CPCTCNT.patientId = CPCTPN.patientId AND CTCT2YN.patientId = CPCTCNT.patientId
  AND CPCTCNT.item = 'FLD.REG.CPCTCNT' AND CPCTCNT.itemValue != ''
  AND CPCTPN.item = 'FLD.REG.CPCTPN' AND CPCTPN.itemValue != ''
  AND CTCT2YN.item = 'FLD.REG.CTCT2YN' AND CTCT2YN.itemValue = 'Yes') a
