import { OrganMetadataCollection } from '../organ/metadata.js';
import { getSampleBlockId, getSampleSectionId, ORGAN_MAPPING } from '../xconsortia/metadata.js';

/**
 * @typedef {object} HubmapMetadata
 * @property {string} uuid
 * @property {string} organ
 * @property {string} organ_source
 * @property {string} assay_type
 * @property {string} dataset_id
 * @property {string} donor_sex
 * @property {string} donor_race
 * @property {string} donor_age
 * @property {string} donor_bmi
 * @property {string} block_sample_id
 * @property {string} section_sample_id
 * @property {string} donor_id
 */

const HUBMAP_ENTITY_ENDPOINT = 'https://entity.api.hubmapconsortium.org/entities/';
const HUBMAP_PORTAL_ENDPOINT = 'https://portal.hubmapconsortium.org/browse/dataset/';
export const ID_KEYWORD = 'hubmap_id';
export const METADATA_FIELDS = [
  'uuid',
  'hubmap_id',
  'origin_samples.organ',
  'dataset_info',
  'mapped_consortium',
  'group_name',
  'group_uuid',
  'donor.mapped_metadata.race',
  'donor.mapped_metadata.sex',
  'donor.mapped_metadata.age_value',
  'donor.mapped_metadata.body_mass_index_value',
  'ancestors',
  'donor.uuid',
];

/**
 * Creates a metadata lookup from the raw metadata
 *
 * @param {object} result Raw metadata
 * @param {OrganMetadataCollection} organMetadata Organ metadata
 */
export function metadataToLookup(result, organMetadata) {
  /** @type {Map<string, HubmapMetadata>} */
  const lookup = new Map();
  for (const hit of result.hits.hits) {
    const {
      _source: {
        hubmap_id,
        uuid,
        origin_samples: [{ organ }],
        dataset_info,
        mapped_consortium,
        group_name,
        group_uuid,
        donor: {
          mapped_metadata: {
            age_value: [donor_age] = [],
            race: [donor_race] = [],
            sex: [donor_sex] = [],
            body_mass_index_value: [donor_bmi] = [],
          } = {},
          uuid: donor_uuid,
        },
        ancestors,
      },
    } = hit;
    const mapped_organ = organMetadata.resolve(ORGAN_MAPPING[organ.toUpperCase()]?.organ_id ?? '');
    const { block_id, rui_location } = getSampleBlockId(ancestors, HUBMAP_ENTITY_ENDPOINT);
    lookup.set(hubmap_id, {
      organ: mapped_organ,
      organ_source: organ,
      uuid,
      assay_type: dataset_info.split('__')[0],
      dataset_id: `${HUBMAP_ENTITY_ENDPOINT}${uuid}`,
      dataset_link: `${HUBMAP_PORTAL_ENDPOINT}${uuid}`,
      dataset_technology: 'OTHER',
      dataset_info,
      consortium_name: mapped_consortium,
      provider_name: group_name,
      provider_uuid: group_uuid,
      donor_id: `${HUBMAP_ENTITY_ENDPOINT}${donor_uuid}`,
      donor_age: donor_age ?? '',
      donor_sex: donor_sex ?? '',
      donor_bmi: donor_bmi ?? '',
      donor_race: donor_race ?? '',
      organ_id: `http://purl.obolibrary.org/obo/UBERON_${mapped_organ.split(':')[1]}`,
      block_id,
      section_id: getSampleSectionId(ancestors, HUBMAP_ENTITY_ENDPOINT),
      rui_location,
    });
  }

  return lookup;
}
