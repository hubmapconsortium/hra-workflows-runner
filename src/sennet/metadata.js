import { getMetadata, getSampleBlockId, getSampleSectionId, ORGAN_MAPPING } from '../xconsortia/metadata.js';

/**
 * @typedef {object} SennetMetadata
 * @property {string} uuid
 * @property {string} organ
 * @property {string} organ_source
 * @property {string} assay_type
 * @property {string} dataset_id
 * @property {string} group_name
 * @property {string} group_uuid
 * @property {string} donor_sex
 * @property {string} donor_race
 * @property {string} donor_age
 * @property {string} donor_bmi
 * @property {string} block_sample_id
 * @property {string} section_sample_id
 * @property {string} donor_id
 */

const SENNET_ENTITY_ENDPOINT = 'https://entity.api.sennetconsortium.org/entities/';
const SENNET_PORTAL_ENDPOINT = 'https://data.sennetconsortium.org/dataset';
export const ID_KEYWORD = 'sennet_id';
export const METADATA_FIELDS = [
  'uuid',
  'sennet_id',
  'origin_sample.organ',
  'dataset_type',
  'group_name',
  'group_uuid',
  'source.source_mapped_metadata.race',
  'source.source_mapped_metadata.sex',
  'source.source_mapped_metadata.age_value',
  'source.source_mapped_metadata.body_mass_index_value',
  'ancestors',
  'source.uuid',
];

/**
 * Creates a metadata lookup from the raw metadata
 *
 * @param {object} result Raw metadata
 */
export function toLookup(result) {
  /** @type {Map<string, SennetMetadata>} */
  const lookup = new Map();
  for (const hit of result.hits.hits) {
    const {
      _source: {
        sennet_id,
        uuid,
        origin_sample: { organ },
        dataset_type: assay_type,
        dataset_info: i,
        group_name,
        group_uuid,
        source: {
          source_mapped_metadata: {
            age_value: [donor_age] = [''],
            race: [donor_race] = [''],
            sex: [donor_sex] = [''],
            body_mass_index_value: [donor_bmi] = [''],
          } = {},
          uuid: donor_uuid,
        },
        ancestors,
      },
    } = hit;
    const mapped_organ = ORGAN_MAPPING[organ.toUpperCase()]?.organ_id ?? '';

    const { block_id, rui_location } = getSampleBlockId(ancestors, SENNET_ENTITY_ENDPOINT);
    lookup.set(sennet_id, {
      organ: mapped_organ,
      organ_source: organ,
      uuid,
      assay_type,
      dataset_id: `${SENNET_ENTITY_ENDPOINT}${uuid}`,
      dataset_link: `${SENNET_PORTAL_ENDPOINT}?uuid=${uuid}`,
      dataset_technology: 'OTHER',
      consortium_name: 'SenNet',
      provider_name: group_name,
      provider_uuid: group_uuid,
      donor_id: `${SENNET_ENTITY_ENDPOINT}${donor_uuid}`,
      donor_age,
      donor_sex,
      donor_bmi,
      donor_race,
      organ_id: `http://purl.obolibrary.org/obo/UBERON_${mapped_organ.split(':')[1]}`,
      block_id,
      section_id: getSampleSectionId(ancestors, SENNET_ENTITY_ENDPOINT),
      rui_location,
    });
  }
  return lookup;
}
