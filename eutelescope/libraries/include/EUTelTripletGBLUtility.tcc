#ifndef EUTelTripletGBLUtility_tcc
#define EUTelTripletGBLUtility_tcc
namespace eutelescope {

template<typename T>
void EUTelTripletGBLUtility::FindTriplets(std::vector<EUTelTripletGBLUtility::hit> const & hits, T const & triplet_sensor_ids, double trip_res_cut, double slope_cut, std::vector<EUTelTripletGBLUtility::triplet> & found_triplets, bool only_best_triplet, bool upstream) {

  if(triplet_sensor_ids.size() != 3){
    throw std::runtime_error("EUTelTripletGBLUtility::FindTriplets called with an invalid set of triplet_sensor_ids (size should be three entries)");
  }

  auto plane0 = static_cast<unsigned>(triplet_sensor_ids[0]);
  auto plane1 = static_cast<unsigned>(triplet_sensor_ids[1]);
  auto plane2 = static_cast<unsigned>(triplet_sensor_ids[2]);

  // get all hit is plane = plane0
  for( auto& ihit: hits ){
    if( ihit.plane != plane0 ) continue; // First plane

    // get all hit is plane = plane2
    for( auto& jhit: hits ){
      if( jhit.plane != plane2 ) continue; // Last plane

      double sum_res_old = -1.;
      // get all hit is plane = plane1
      for( auto& khit: hits ){
	if( khit.plane != plane1 ) continue; // Middle plane

	// Create new preliminary triplet from the three hits:
	EUTelTripletGBLUtility::triplet new_triplet(ihit,khit,jhit);

    //Create triplet slope plots
    if(upstream == 1){
		upstreamTripletSlopeX->fill(new_triplet.getdx()*1E3/new_triplet.getdz()); //factor 1E3 to convert from rad to mrad. To be checked
		upstreamTripletSlopeY->fill(new_triplet.getdy()*1E3/new_triplet.getdz());
	} else {
		downstreamTripletSlopeX->fill(new_triplet.getdx()*1E3/new_triplet.getdz());
		downstreamTripletSlopeY->fill(new_triplet.getdy()*1E3/new_triplet.getdz());
	}
	// Setting cuts on the triplet track angle:
	if( fabs(new_triplet.getdx()) > slope_cut * new_triplet.getdz()) continue;
	if( fabs(new_triplet.getdy()) > slope_cut * new_triplet.getdz()) continue;
    
    //Create triplet residual plots
    if(upstream == 1){
		upstreamTripletResidualX->fill(new_triplet.getdx(plane1));
		upstreamTripletResidualY->fill(new_triplet.getdy(plane1));
	} else {
		downstreamTripletResidualX->fill(new_triplet.getdx(plane1));
		downstreamTripletResidualY->fill(new_triplet.getdy(plane1));
	}
	// Setting cuts on the triplet residual on the middle plane
	if( fabs(new_triplet.getdx(plane1)) > trip_res_cut) continue;
	if( fabs(new_triplet.getdy(plane1)) > trip_res_cut) continue;

    // Edo: This (best triplets) is hardcoded as false in EUTelAlignGBL. Is this really useful?
    if(only_best_triplet) {
		// For low threshold (high noise) and/or high occupancy, use only the triplet with the smallest sum of residuals on plane1
		double sum_res = sqrt(new_triplet.getdx(plane1)*new_triplet.getdx(plane1) + new_triplet.getdy(plane1)*new_triplet.getdy(plane1));
		if(sum_res < sum_res_old){
	
		  // Remove the last one since it fits worse, not if its the first
		  found_triplets.pop_back();
		  // The triplet is accepted, push it back:
		  found_triplets.emplace_back(new_triplet);
		  streamlog_out(DEBUG2) << new_triplet;
		  sum_res_old = sum_res;
		}

		// update sum_res_old on first iteration
		if(sum_res_old < 0.) {
		  // The triplet is accepted, push it back:
		  found_triplets.emplace_back(new_triplet);
		  streamlog_out(DEBUG2) << new_triplet;
		  sum_res_old = sum_res;
		}
	} else {	
		found_triplets.emplace_back(new_triplet);
	}
      }//loop over hits
    }//loop over hits
  }// loop over hits

  //return triplets;
}

/********* 
    * current_multiplets : current set of the multiplets built with the already used planes
    * in_hits : all hits in the telescope for this event. The two indices are [ID plane][ID hits on the plane]
    ******/
std::vector<eutelescope::EUTelTripletGBLUtility::multiplet> EUTelTripletGBLUtility::RecursiveMultipletBuilding(std::vector<eutelescope::EUTelTripletGBLUtility::multiplet> &current_multiplets, std::vector< std::vector<eutelescope::EUTelTripletGBLUtility::hit> >& in_hits)
{

    if(in_hits.size() == 2) { //Initial condition
        std::vector<EUTelTripletGBLUtility::hit> all_hits_p0 = in_hits.at(0);
        std::vector<EUTelTripletGBLUtility::hit> all_hits_p1 = in_hits.at(1);
        std::vector<EUTelTripletGBLUtility::multiplet> all_output_multiplets;
        
        for(EUTelTripletGBLUtility::hit hit_p0 : in_hits.at(0)) 
        for(EUTelTripletGBLUtility::hit hit_p1 : in_hits.at(1))
        {
            EUTelTripletGBLUtility::multiplet out_multiplet;
            std::vector<EUTelTripletGBLUtility::hit> current_multiplet_hits;
            current_multiplet_hits.push_back(hit_p0);
            current_multiplet_hits.push_back(hit_p1);
            //--- TODO: Add some cuts here !
            out_multiplet.fillmultiplet(current_multiplet_hits);
            all_output_multiplets.push_back(current_multiplet_hits);
        }
        return all_output_multiplets;
    } else { //Recursivityzation
        std::vector<EUTelTripletGBLUtility::hit> current_plane_hits = in_hits.at(1);
        in_hits.erase(in_hits.begin() + 1);
            
        std::vector<EUTelTripletGBLUtility::multiplet> sub_multiplet = RecursiveMultipletBuilding(current_multiplets, in_hits);
        
        for(EUTelTripletGBLUtility::hit current_hit : current_plane_hits)
        for(EUTelTripletGBLUtility::multiplet current_build_multiplet : sub_multiplet)
        {
            current_build_multiplet.add_hit(current_hit);
            //--- TODO : Add some cuts here
            //--- TODO Maybe : Add an option to check whether a plane is DUT ?
        }
    }
}

template<typename T>
void EUTelTripletGBLUtility::FindMultiplets(std::vector<EUTelTripletGBLUtility::hit> const & hits, std::vector<T> const & multiplet_sensor_ids, double multip_res_cut, double multip_slope_cut, std::vector<EUTelTripletGBLUtility::multiplet> & found_multip, bool only_best_multiplet) {

  vector<unsigned> planes;
  for(unsigned current = 0 ; current < multiplet_sensor_ids.size() ; current++) 
    planes.push_back( static_cast<unsigned> (multiplet_sensor_ids[current]) );
    
  std::vector< std::vector<EUTelTripletGBLUtility::hit> > hits_by_plane;
  for(auto current_plane : planes) 
    for(auto current_hit : hits) 
      if(current_hit.plane == current_plane) hits_by_plane.push_back(current_hit);

  // get all hit is plane = plane0
  for( auto& ihit: hits ){
    if( ihit.plane != plane0 ) continue; // First plane

    // get all hit is plane = plane2
    for( auto& jhit: hits ){
      if( jhit.plane != plane2 ) continue; // Last plane

      double sum_res_old = -1.;
      // get all hit is plane = plane1
      for( auto& khit: hits ){
	if( khit.plane != plane1 ) continue; // Middle plane

	// Create new preliminary triplet from the three hits:
	EUTelTripletGBLUtility::triplet new_triplet(ihit,khit,jhit);

    //Create triplet slope plots
    if(upstream == 1){
		upstreamTripletSlopeX->fill(new_triplet.getdx()*1E3/new_triplet.getdz()); //factor 1E3 to convert from rad to mrad. To be checked
		upstreamTripletSlopeY->fill(new_triplet.getdy()*1E3/new_triplet.getdz());
	} else {
		downstreamTripletSlopeX->fill(new_triplet.getdx()*1E3/new_triplet.getdz());
		downstreamTripletSlopeY->fill(new_triplet.getdy()*1E3/new_triplet.getdz());
	}
	// Setting cuts on the triplet track angle:
	if( fabs(new_triplet.getdx()) > slope_cut * new_triplet.getdz()) continue;
	if( fabs(new_triplet.getdy()) > slope_cut * new_triplet.getdz()) continue;
    
    //Create triplet residual plots
    if(upstream == 1){
		upstreamTripletResidualX->fill(new_triplet.getdx(plane1));
		upstreamTripletResidualY->fill(new_triplet.getdy(plane1));
	} else {
		downstreamTripletResidualX->fill(new_triplet.getdx(plane1));
		downstreamTripletResidualY->fill(new_triplet.getdy(plane1));
	}
	// Setting cuts on the triplet residual on the middle plane
	if( fabs(new_triplet.getdx(plane1)) > trip_res_cut) continue;
	if( fabs(new_triplet.getdy(plane1)) > trip_res_cut) continue;

    // Edo: This (best triplets) is hardcoded as false in EUTelAlignGBL. Is this really useful?
    if(only_best_triplet) {
		// For low threshold (high noise) and/or high occupancy, use only the triplet with the smallest sum of residuals on plane1
		double sum_res = sqrt(new_triplet.getdx(plane1)*new_triplet.getdx(plane1) + new_triplet.getdy(plane1)*new_triplet.getdy(plane1));
		if(sum_res < sum_res_old){
	
		  // Remove the last one since it fits worse, not if its the first
		  found_triplets.pop_back();
		  // The triplet is accepted, push it back:
		  found_triplets.emplace_back(new_triplet);
		  streamlog_out(DEBUG2) << new_triplet;
		  sum_res_old = sum_res;
		}

		// update sum_res_old on first iteration
		if(sum_res_old < 0.) {
		  // The triplet is accepted, push it back:
		  found_triplets.emplace_back(new_triplet);
		  streamlog_out(DEBUG2) << new_triplet;
		  sum_res_old = sum_res;
		}
	} else {	
		found_triplets.emplace_back(new_triplet);
	}
      }//loop over hits
    }//loop over hits
  }// loop over hits

  //return triplets;
}

}//namespace
#endif
