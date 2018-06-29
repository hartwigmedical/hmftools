hmf_id_hashes has been generated using v1.0 of hmf-id-generator.

The set of 3318 patient ids was queried on June 29th 2018, 10:21:00 with the following:

select distinct patientIdentifier from
(select patientIdentifier from patient
UNION
SELECT concat('CPCT02', CPCTCNT.itemValue, LPAD(RIGHT(CPCTPN.itemValue,4), 4, '0')) as patientIdentifier
 FROM drupEcrf CTCT2YN,  drupEcrf CPCTCNT, drupEcrf CPCTPN
WHERE CPCTCNT.patientId = CPCTPN.patientId AND CTCT2YN.patientId = CPCTCNT.patientId
  AND CPCTCNT.item = 'FLD.REG.CPCTCNT' AND CPCTCNT.itemValue != ''
  AND CPCTPN.item = 'FLD.REG.CPCTPN' AND CPCTPN.itemValue != ''
  AND CTCT2YN.item = 'FLD.REG.CTCT2YN' AND CTCT2YN.itemValue = 'Yes') a

