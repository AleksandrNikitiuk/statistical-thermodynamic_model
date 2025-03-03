function [energy_barrier] = get_energy_barrier(free_energy)

local_minima = free_energy(islocalmin(free_energy));

if isempty(local_minima)
  local_maxima = free_energy(islocalmax(free_energy));
  if isempty(local_maxima)
    [~,n_maximum] = max(free_energy);
    if n_maximum == 1
      energy_barrier = 0;
    else
      energy_barrier = free_energy(n_maximum) - min(free_energy);
    end
  else
    energy_barrier = local_maxima(1) - min(free_energy(free_energy < local_maxima(1)));
  end
else
  local_maxima = free_energy(islocalmax(free_energy));
  if isempty(local_maxima)
    free_energy = free_energy(find(free_energy == local_minima(1)):end);
    [~,n_maximum] = max(free_energy);
    if free_energy(n_maximum) == local_minima(1)
      energy_barrier = 0;
    else
      energy_barrier = free_energy(n_maximum) - local_minima(1);
    end
  else
    energy_barrier = local_maxima(1) - local_minima(1);
  end
end

end