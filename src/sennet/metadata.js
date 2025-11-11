import { OrganMetadataCollection } from '../organ/metadata.js';
import { getSampleBlockId, getSampleSectionId, ORGAN_MAPPING } from '../xconsortia/metadata.js';

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
const SENNET_ANCESTORS_ENDPOINT = 'https://entity.api.sennetconsortium.org/ancestors/';
const SENNET_PORTAL_ENDPOINT = 'https://data.sennetconsortium.org/dataset';
export const ID_KEYWORD = 'sennet_id';
export const METADATA_FIELDS = [
  'uuid',
  'sennet_id',
  'origin_samples.organ',
  'dataset_info',
  'group_name',
  'group_uuid',
  'sources.uuid',
];

/**
 * Creates a metadata lookup from the raw metadata
 *
 * @param {object} result Raw metadata
 * @param {OrganMetadataCollection} organMetadata Organ metadata
 */
export async function toLookup(result, organMetadata, token = undefined) {
  // Headers to use in authenticated fetch requests
  const headers = {};
  if (token) {
    headers['Authorization'] = `Bearer ${token}`;
  }

  /** @type {Map<string, SennetMetadata>} */
  const lookup = new Map();
  for (const hit of result.hits.hits) {
    const {
      _source: {
        sennet_id,
        uuid,
        origin_samples: [{ organ }],
        dataset_info,
        group_name,
        group_uuid,
        sources: [{ uuid: donor_uuid }],
      },
    } = hit;
    const source = await fetch(`${SENNET_ENTITY_ENDPOINT}${donor_uuid}`, { headers }).then((r) => r.json());
    const {
      source_mapped_metadata: {
        age: { value: [donor_age] = [''] } = {},
        race: { value: [donor_race] = [''] } = {},
        sex: { value: [donor_sex] = [''] } = {},
        body_mass_index: { value: [donor_bmi] = [''] } = {},
      } = {},
    } = source;

    const ancestors = await fetch(`${SENNET_ANCESTORS_ENDPOINT}${uuid}`, { headers }).then((r) => r.json());
    if (ancestors.error) {
      console.error(`Error getting ancestors for ${uuid}: ${ancestors.error}`);
      throw new Error(ancestors.error);
    }

    const mapped_organ = organMetadata.resolve(ORGAN_MAPPING[organ.toUpperCase()]?.organ_id ?? organ);
    const { block_id, rui_location } = getSampleBlockId(ancestors, SENNET_ENTITY_ENDPOINT);
    lookup.set(sennet_id, {
      organ: mapped_organ,
      organ_source: organ,
      uuid,
      assay_type: dataset_info.split('__')[0],
      dataset_id: `${SENNET_ENTITY_ENDPOINT}${uuid}`,
      dataset_link: `${SENNET_PORTAL_ENDPOINT}?uuid=${uuid}`,
      dataset_technology: 'OTHER',
      dataset_info,
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
